function [mcs, status, obj] = CNAmultiImposedConstrISfinder(cnap,T,t,D,d,koCost,kiCost,maxNumInterv,timeLimit,maxSolutions,use_compression,use_bigM,enum_method,verbose)
% function 'CNAmultiImposedConstrISfinder'
% ---------------------------
% Given a mass-flow project and a set of 'undesired'
% (target) flux vectors (defined by one or multiple matrices T and vectors t) and (optionally) a set of
% 'desired' flux vectors (defined by one or multiple matrices D and vectors d) an intervention strategy is
% computed fulfilling the following properties: Knocking out the reactions of the MCS
% (i.e., setting the corresponding rate to zero) ensures that all target flux vectors
% v obeying
%
%       cnap.stoichiMat * v = 0
%       reac.Min <= v <= v.reacMax
%       T*v <= t
%
% will be blocked (are infeasible) whereas at least one flux vector r fulfilling
%
%       cnap.stoichiMat * r = 0
%       cnap.reacMin <= r <= cnap.reacMax
%       D*r <= d
%
% will be kept functional, possibly by activating other reactions
%    (if D is empty, no such flux vector must exist).
%
% Usage:
%
% Input:
%   cnap: a CellNetAnalyzer project
%   T: A cell array of matrices specifying (with vector t) the target flux
%      regions as given above.  T has Dimension numTargetConst(trgt_idx) x cnap.numr
%
%   t: A cell array of vectors specifying (with matrix T) the target flux
%      regions as given above.  t has Dimension numTargetConst(trgt_idx) x 1
%
%   D: One matrix or a cell array of matrices specifying (with vector d) the desired flux
%      regions as given above.  D has Dimension numDesiredConst(des_idx) x cnap.numr
%
%   d: One vector or a cell array of vectors specifying (with matrix D) the desired flux
%      regions as given above.  d has Dimension numDesiredConst(des_idx) x 1
%
%   notknockable: row vector with indices of reactions which cannot be cut (knocked-out).
%
%   max_num_interv: (optional; default: []) maximum size (cardinality) the IS
%      may have, use [] if you do not want to limit IS size
%
%   time_limit: maximal time after which to stop calculation
%
%   max_solutions: maximum number of solutions that are found
%
% Output:
%   is: intervention strategies
%
if nargin < 11 || isempty(use_compression)
    use_compression = 1;
end
if nargin < 12 || isempty(use_bigM)
    use_bigM = 0;
end
if nargin < 13 || isempty(enum_method)
    enum_method = 1;
end
if nargin < 14 || isempty(verbose)
    verbose = 1;
end
if isempty(kiCost)
    kiCost = nan(cnap.numr,1);
end
if isempty(D)
    lb_D={};
    ub_D={};
end
disp('== MCS Computation ==');
c_macro = cnap.macroDefault;


% Transform ki & ko vector to row vector with 'nan' as notknockable and a double as 'weight'
kiCost = kiCost(:)';
koCost = koCost(:)';
koCost(~isnan(kiCost)) = nan; % knock-ins 'override' knock-outs

% Bound the system, otherwise MCS can be found that disrupt desired region/s
if ~isempty(D)
    cnap.reacMax(cnap.reacMax ==  inf) =  1e5;
    cnap.reacMin(cnap.reacMin == -inf) = -1e5;
end

% parallel or not depends on the problem size
if cnap.numr > 500 ,parforArg = inf; else, parforArg = 0; end

%% test feasibility of target and desired vectors
displ('Verifying that D and T regions are feasible in the original model.',verbose);
cnap.objFunc = zeros(cnap.numr,1);
testRegionFeas(cnap,c_macro,T,t,D,d);

%% Determine boundaries to whole model
displ('FVA to determine model bounds.',verbose);
[cnap.reacMin, cnap.reacMax] = CNAfluxVariability(cnap,[],c_macro,2);

%% Determine boundaries to desired scenarios and identify further notknockables
displ('FVA to determine model bounds under desired constraints.',verbose);
for i = 1:length(D)
    [lb_D{i}, ub_D{i}] = CNAfluxVariability(cnap,[],c_macro,2,1:cnap.numr,D{i},d{i});
    lb_D{i}( abs(lb_D{i})<=cnap.epsilon ) = 0;
    ub_D{i}( abs(ub_D{i})<=cnap.epsilon ) = 0;
    essential = (sign(lb_D{i}).*sign(ub_D{i})) == 1;
    koCost(essential) = nan; % make essential reactions "notknockable"
end
if ~isempty(D) && any(any(isinf([cell2mat(lb_D) cell2mat(ub_D)]),2))
    error(['Model must be bound to use desired constraints. But the follwing reactions are unbound: '...
        strjoin(cellstr(cnap.reacID(any(isinf([cell2mat(lb_D) cell2mat(ub_D)]),2),:)),', ')])
end

%% compress network if indicated
if use_compression
    kiCost_full = kiCost;
    koCost_full = koCost;
    maxNumInterv_full = maxNumInterv;
    %% replace the full model variables by the compressed ones
    displ('Compressing model.',verbose);
    [cnap, T, D, koCost, kiCost,  cmp_transf_mat] = ...
        compress(cnap,T,D,koCost,kiCost);
    for i = 1:length(D)
        [lb_D{i}, ub_D{i}] = generateReducedBounds(lb_D{i},ub_D{i},cmp_transf_mat);
    end
    %% test again feasibility of target and desired vectors
    displ('Verifying D and T region feasibility in compressed model.',verbose);
    testRegionFeas(cnap,c_macro,T,t,D,d);
else
    cmp_transf_mat = eye(cnap.numr);
end

%% bigM
if use_bigM
    M = max(max(abs([cell2mat(cnap.reacMin) cell2mat(cnap.reacMax)]))) * 1.1; % Choose bigM 10% bigger than largest bound
else
    M = 0;
end
%% Make sure knock-in-able reactions have boundaries that include 0
for kis = find(~isnan(kiCost))
    for i = 1:length(D)
        if     lb_D{i}(kis) < 0 && ub_D{i}(kis) < 0
            ub_D{i}(kis) = 0;
        elseif lb_D{i}(kis) > 0 && ub_D{i}(kis) > 0
            lb_D{i}(kis) = 0;
        end
    end
    if     cnap.reacMin(kis) < 0 && cnap.reacMax(kis) < 0
        cnap.reacMax(kis) = 0;
    elseif cnap.reacMin(kis) > 0 && cnap.reacMax(kis) > 0
        cnap.reacMin(kis) = 0;
    end
end

%% integrate flux bounds in Target system
%P Axel: Maybe replace positive and negative bounds by +/-inf
for i = 1:length(T)
    Tidx=size(T{i},1);
    for h = 1:cnap.numr
        if cnap.reacMin(h)~=0 && ~isinf(cnap.reacMin(h))
            Tidx=Tidx+1;
            T{i}(Tidx,h)=-1;
            t{i}(Tidx)=-cnap.reacMin(h);
        end
        if ~isinf(cnap.reacMax(h))
            Tidx=Tidx+1;
            T{i}(Tidx,h)=1;
            t{i}(Tidx)=cnap.reacMax(h);
        end
    end
end

%% compute intervention strategies
N  = initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
obj= ConstrainedMinimalCutSetsEnumerator(  N, ...
    cnap.reacMin' >= 0, [], ...  Irr
    T, t, ...                 Target
    {~isnan(koCost), full(sparse(1,find(~isnan(koCost)),koCost(~isnan(koCost)),1,length(koCost)))}, ... KO
    lb_D, ...        Bounds
    ub_D, ...
    D, d, ...        Desired
    {~isnan(kiCost), full(sparse(1,find(~isnan(kiCost)),kiCost(~isnan(kiCost)),1,length(kiCost)))},...  KI
    M);
cplex_inner= setup_cplex_inner_class_access();
if isunix
    slurm_job_mem = str2double(getenv('SLURM_REQMEM'));
    if ~isnan(slurm_job_mem)
        usemem = round(slurm_job_mem*0.85); % use 85% of allocated memory
        delete(gcp('nocreate')); % stop parpool to free additional memory
    else
        [~,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        usemem = round(stats(6)/1e3*0.6); % use 60% of available memory
    end
else
    [~,mem] = memory;
    usemem = mem.PhysicalMemory.Available/1e6*0.9; % use up to 90% of currently available memory
end
cpxparam = obj.cpx.getParameterSet;
cpxparam.setParam(cplex_inner.DoubleParam.WorkMem, usemem); % Allocate memory
if ~ispc
cpxparam.setParam(cplex_inner.IntParam.RootAlg, 2); % Use Solution method that doesn't lead to crash
end
cpxparam.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Feasibility);
cpxparam.setParam(cplex_inner.IntParam.FPHeur,1);
% cpxparam.setParam(cplex_inner.IntParam.Probe,3);
cpxparam.setParam(cplex_inner.IntParam.MIPDisplay, 2); % Ignored when setOut = []
%cpxparam.setParam(cplex_inner.IntParam.NodeSel, 0); % Use "Dept first" node selection.
                                    % Leads to longer runtimes, but consumes much less memory.
                                    % When first integer solution is needed, this is an okay
                                    % tradeoff.

slurm_job_cpus = str2double(getenv('SLURM_CPUS_ON_NODE'));
if ~isnan(slurm_job_cpus)
    cpxparam.setParam(cplex_inner.IntParam.Threads,slurm_job_cpus); % Use all allocated SLURM cores.
    % cpxparam.setParam(cplex_inner.IntParam.Threads,min(6,slurm_job_cpus)); % Use all allocated SLURM cores, but a maximum of 6 cores.
else
    cpxparam.setParam(cplex_inner.IntParam.Threads,min(4,(2*feature('numcores'))-1)); % Use one core less than available, but a maximum of 4 cores.
end
obj.cpx.setParameterSet(cpxparam);
obj.cpx.setOut([]);

knockable = arrayfun(@(x) x.getUB(),obj.z_vars);
text = ['Starting MCS computation. Problem size (columns): ' num2str(obj.cpx.getNcols) '. Knockable reactions: ' num2str(sum(knockable~=0))];
if ~isinf(timeLimit),  text = [text '. Time limit: ' num2str(timeLimit)];  end
displ(text,verbose);
displ(['Using up to ' num2str(usemem) ' MB of memory and ' num2str(cpxparam.getParam(cplex_inner.IntParam.Threads)) ' threads'],verbose);
starttime = now;
if enum_method == 1
    status = 0;
    mcs = ones(cnap.numr,0);
    ivCost = kiCost;
    ivCost(~isnan(koCost)) = koCost(~isnan(koCost));
    ivCost(isnan(ivCost))  = 0;
    obj.evs_sz.setUB(maxNumInterv);
    endtime = now*86400 + timeLimit;
    %% compute smallest mcs
    if sum(knockable) <= 100 % in larger networks
        disp('Computing smallest MCS to find the lower bound for the MCS size.');
        [mcs, status] = findMCS(obj, 'minimize', cplex_inner, endtime-now*86400);
        if status == 0
            obj.evs_sz.setLB(koCost(~isnan(koCost))*mcs(~isnan(koCost)));
            obj= add_exclusion_constraints(obj, mcs, sum(mcs));
        end
    end
    %% compute largest cut set (only for complete enumeration)
    if sum(knockable) <= maxNumInterv && maxSolutions == inf && status == 0
        disp('Computing largest cut set to define an upper bound for the MCS size.');
        [sol, status] = findMCS(obj, 'maximize', cplex_inner, endtime-now*86400);
        if status == 0
            maxNumInterv = max(sum(mcs.*ivCost'));
            arrayfun(@(x,y) x.setUB(y),obj.z_vars,sol);
            while status == 0 && size(mcs,2) < maxSolutions
                [sol, status] = findMCS(obj, 'minimize', cplex_inner, endtime-now*86400);
                if  status == 0
                    mcs = [mcs, sol];
                    obj= add_exclusion_constraints(obj, sol, sum(sol));
                end
            end
            arrayfun(@(x,y) x.setUB(y),obj.z_vars,knockable);
            status = 0;
        end
    end
    %% search for MCS iteratively
    disp('Computing MCS (iteratively).');
    while size(mcs,2) < maxSolutions && status == 0
        [sol, status] = findMCS(obj, 0, cplex_inner, endtime-now*86400);
        if status == 0
            arrayfun(@(x,y) x.setUB(y),obj.z_vars,sol);
            while status == 0 && size(mcs,2) < maxSolutions
                [sol, status] = findMCS(obj, 'minimize', cplex_inner, endtime-now*86400);
                if  status == 0 && any(sol)
                    mcs = [mcs, sol];
                    fprintf('.');
                    obj= add_exclusion_constraints(obj, sol, sum(sol));
                    if verbose % output
                        text = [num2str(round(86400*(now - starttime)*100)/100,'%.2f') ' seconds : '];
                        weighted_mcs = sum(mcs.*repmat(ivCost',1,size(mcs,2)),1);
                        for mcs_size = unique(sum(weighted_mcs,1))
                            text = [text, num2str(sum(weighted_mcs==mcs_size)), ' of size ', num2str(mcs_size) ,';  '];
                        end
                        disp(text);
                    end
                end
            end
            arrayfun(@(x,y) x.setUB(y),obj.z_vars,knockable);
            status = 0;
        end
    end
    fprintf(char(10));
elseif enum_method==2
    [obj, mcs, status]= shortest_minimal_cut_sets(obj, maxNumInterv, maxSolutions, timeLimit, 2);
    if status.equals(cplex_inner.Status.Feasible) || status.equals(cplex_inner.Status.Optimal)
        status = 0;
    elseif ~isinf(timeLim) && obj.getCplexStatus().equals(cplex_inner.CplexStatus.AbortTimeLim)
        status = 1;
    elseif status.equals(cplex_inner.Status.Infeasible)
        status = 2;
    end
end
if size(mcs,2) == maxSolutions
    disp('Solution pool limit reached.');
end

if ~isnan(mcs)
    %% decompress intervention strategies
    if use_compression
        compSols = size(mcs,2);
        displ('Decompress MCS.',verbose);
        kiCost = kiCost_full;
        koCost = koCost_full;
        ivCost = kiCost;
        ivCost(~isnan(koCost)) = koCost(~isnan(koCost));
        if maxNumInterv_full == inf
            maxNumInterv_full = 1e9; % workaround
        end
        ivCost(isnan(ivCost))  = maxNumInterv_full+1; % higher number, so that mcs containing this intervention are sortet out
        mcs = expand_mcs(mcs, cmp_transf_mat');
        % check if knock-out and knock-in costs are still met
        mcs_affordable = sum(mcs.*repmat(ivCost',1,size(mcs,2)),1) <= maxNumInterv_full;
        mcs = mcs(:,mcs_affordable);
    end

    %% set knock-in/out canditates to [-1: knocked out] or [1: knocked in] or [nan: not knocked in]
    mcs = -double(mcs);
    mcs(~isnan(kiCost),:) = double(-mcs(~isnan(kiCost),:));
    for i = find(~isnan(kiCost(:)))'
        mcs(i, mcs(i,:)==0) = nan;
    end
    
    if size(mcs,2) > 0 && status == 2
        status = 0; 
    elseif status == 1
        disp('Time limit reached.');
    elseif status == 0

    displ(['Total computation time: ' num2str(86400*(now - starttime))],verbose);
    if use_compression, text =  ['. (' num2str(compSols) ' compressed)']; else , text = []; end
        if all(all(isnan(mcs)))
            displ('no MCS found',verbose);
        else
            displ(['MCS found: ' num2str(size(mcs,2)) text],verbose);
        end
    end
else
    disp('Problem was infeasible.');
end
end

%% Functions:
%% 1. compress
function [cmp_cnap, cmp_T, cmp_D, cmp_koCost, cmp_kiCost, cmp_mapReac] = compress(cnap,T,D,koCost,kiCost)
    %% Prepare
    % identify essential reactions and adapt notknockable vector
    non_compress_reacs = any([cell2mat(D') ; cell2mat(T')],1);
    %         non_compress_reacs( isnan(koCost)) = true; % don't compress notknockable %P Can be compressed, now
    non_compress_reacs(~isnan(kiCost)) = true; % don't compress knockinable
    non_compress_reacs = find(non_compress_reacs);

    r_off      = abs(cnap.reacMin)<=cnap.epsilon & ...
        abs(cnap.reacMax)<=cnap.epsilon;

    javastderr= java.lang.System.err;
    java.lang.System.setErr(java.io.PrintStream('cplex_stderr.log'));
    %% Compress
    [dum1,dum2,cmp_mapReac,dum3,cmp_cnap]=CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,r_off,1);
    java.lang.System.setErr(javastderr);
    % remap IS-search parameters
    cmp_T = cellfun(@(x) x*cmp_mapReac,T,'UniformOutput',0);
    cmp_D = cellfun(@(x) x*cmp_mapReac,D,'UniformOutput',0);
    lumpedReacs = double(cmp_mapReac ~= 0);
    lumpedReacs(lumpedReacs == 0) = nan;
    cmp_koCost = min(lumpedReacs.*repmat(koCost',1,size(lumpedReacs,2)));
    cmp_kiCost = arrayfun(@(x) sum(kiCost(~isnan(lumpedReacs(:,x)))),1:size(lumpedReacs,2)); % for KI, all lumped reactions need to be knocked in
end

%% 2. generate reduced bounds
function [cmp_lb, cmp_ub] = generateReducedBounds(lb,ub,map)
    [cmp_lb, cmp_ub] = deal(nan(size(map,2),1));
    for k=1:size(map,2)
        forw=find(map(:,k)>0);
        revs=find(map(:,k)<0);
        if(~isempty(forw))
            cmp_lb(k) = max(lb(forw)./map(forw,k));
            cmp_ub(k) = min(ub(forw)./map(forw,k));
        end
        if(~isempty(revs))
            cmp_lb(k) = max([cmp_lb(k);ub(revs)./map(revs,k)]);
            cmp_ub(k) = min([cmp_ub(k);lb(revs)./map(revs,k)]);
        end
    end
end
%% 3. optimization
function [solution, status] = findMCS(obj, sense, cp_inner, timeLim)
    solution = nan;
    if timeLim<=0
        status = 1;
        return;
    end
    if ~isempty(obj.cpx.getObjective()) , obj.cpx.remove(obj.cpx.getObjective()); end
    if ~isinf(timeLim), obj.cpx.setParam(cp_inner.DoubleParam.TiLim,timeLim);
    end
    for j = 0:(obj.cpx.getNMIPStarts-1), obj.cpx.deleteMIPStarts(0); % delete MIP starts
    end
    switch sense
        case 'maximize'
            obj.cpx.addMaximize().setExpr(obj.obj_expr);
        case 'minimize'
            obj.cpx.addMinimize().setExpr(obj.obj_expr);
        otherwise
    end
    obj.solve();
    status = obj.getStatus();
    if status.equals(cp_inner.Status.Feasible) || status.equals(cp_inner.Status.Optimal)
        solution = retrieve_optimal_solution(obj);
        status = 0;
    elseif ~isinf(timeLim) && obj.getCplexStatus().equals(cp_inner.CplexStatus.AbortTimeLim)
        status = 1;
    elseif status.equals(cp_inner.Status.Infeasible)
        status = 2;
    end
    if ~isempty(obj.cpx.getObjective()) , obj.cpx.remove(obj.cpx.getObjective()); end
end
%% 4. adding exclusion constraints
function obj= add_exclusion_constraints(obj, zv, support_size)
    if obj.split_z
        zv= [zv; zv];
    end
    if support_size == 0
        disp('Empty exclusion constraint detected.');
        costs = [];
        iter = obj.evs_sz.getExpr.linearIterator;
        while iter.hasNext
            iter.next;
            costs = [costs iter.getValue];
        end
        obj.evs_sz.setLB(min(costs));
    else
        m= obj.mipmat.addRows(zeros(1, size(zv, 2)), support_size - 1, [], []);
        [k, j, val]= find(zv);
        obj.mipmat.setNZs(m + j - 1, obj.first_z + k - 1, val); %A Java indices, implicit transposition
    end
end
%% 5. test region feasibility
function testRegionFeas(cnap,c_macro,T,t,D,d)
    for i = 1:length(t)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, T{i}, t{i}))
            error(['At least one target region (T' num2str(i) ') is infeasible in the original model']);
        end
    end
    for i = 1:length(d)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, D{i}, d{i}))
            error(['At least one desired region (D' num2str(i) ') is infeasible in the original model']);
        end
    end
end
% Suppress code warnings
function displ(x,flag)
    if flag, disp(x); end
end
%#ok<*ASGLU>
%#ok<*AGROW>