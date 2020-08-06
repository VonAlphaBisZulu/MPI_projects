function [mcs,reacNames,reg, obj, cs]=CNAfindRegMCS2(cnap,T,t,D,d,notknockable,maxMCSnumPre,maxMCSnum,maxMCSsize,filename,reac_off,regulation,time_limit1,time_limit2,EpGap,inverseCutCount)
%
% Similar to 'CNAregMCSEnumerator', but only attempts to find a single (not necessarily smallest) cMCS.
%
% -->  Computes Constrained Minimal Cut Sets (cMCSs) in mass-flow networks.
%      Optionally, up- and down-regulation of certain reaction fluxes may be
%      considered as intervention strategies in combination with
%      reaction cuts (knockouts) yielding then 'regulatory cMCSs'.
%
%  Usage: [mcs,reacNames,reg,sys, obj, cs]=CNAfindRegMCS(cnap,T,t,D,d,notknockable,
%						maxMCSsize,filename,reac_off,regulation,time_limit,EpGap)
%
% Given a mass-flow project and a set of 'undesired'
% (target) flux vectors (defined by matrix T and vector t) and (optionally) a set of
% 'desired' flux vectors (defined by matrix D and vector d) a Minimal Cut Set is
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
% will be kept functional (if D is empty, no such flux vector must exist).
%
% In order to run this function it is necessary that both the MATLAB CPLEX and Java CPLEX
% interfaces work. If CPLEX is installed under /cluster/apps/cplex-124 the commands for
% this are:
%      addpath('/cluster/apps/cplex-124/cplex/matlab/');
%      javaaddpath('/cluster/apps/cplex-124/cplex/lib/cplex.jar');
% Additionally, the MATLAB JVM needs to have the CPLEX shared library on its library path i
% which must be set up before (!!) starting MATLAB. For MATLAB versions up to 2010b
% this can be achieved by adding
%       /cluster/apps/cplex-124/cplex/bin/x86-64_sles10_4.1
% to Matlab's librarypath.txt (or javalibrarypath.txt) configuration file (see also manual).
%
%
% Inputs (the first 5 arguments are mandatory):
% ---------------------------------------------
%
%   cnap: (mandatory) is a CellNetAnalyzer (mass-flow) project variable representing
%         a (metabolic) reaction network. You may easily generate such a structure
%         by the CNAgenerateMFNetwork function.
%
%         The function accesses the following fields of cnap (see manual):
%
%      cnap.stoichmat: the stoichiometric matrix of the network
%      cnap.numr: number of reactions (columns in cnap.stoichMat)
%      cnap.numis: number of internal species
%      cnap.mue: index of the biosynthesis reaction; can be empty
%      cnap.macroComposition: matrix defining the stoichiometry of
%        the macromolecules with respect to the metabolites (species); matrix
%        element macroComposition(i,j) stores how much of metabolite i (in mmol)
%        is required to synthesize 1 gram of macromolecule j
%      cnap.specInternal: vector with the indices of the internal species
%      cnap.macroDefault: default concentrations of the macromolecules
%      cnap.reacMin: lower boundaries of reaction rates
%       (if reacMin(i)=0 --> reaction i is irreversible)
%      cnap.reacMax: upper boundaries of reaction rates
%
%   T: the matrix specifying (with vector t) the target flux vectors as given above.
%      T has Dimension numTargetConst x cnap.numr
%
%   t: the vector specifying (with matrix T) the target flux vectors as given above.
%      t has Dimension numTargetConst x 1
%
%   D: the matrix specifying (with vector d) the desired flux vectors as given above.
%      D has Dimension numDesiredConst x cnap.numr. D and d can be empty
%
%   d: the vector specifying (with matrix D) the desired flux vectors as given above.
%      d has Dimension numDesiredConst x 1. D and d can be empty
%
%   notknockable: row vector with indices of reactions which cannot be cut (knocked-out).
%
%   maxMCSnumPre: (optional; default: []) maximum size (cardinality) the cMCS
%      may have, use [] if you do not want to limit cMCS size
%
%   maxMCSnum: (optional; default: []) maximum size (cardinality) the cMCS
%      may have, use [] if you do not want to limit cMCS size
%
%   maxMCSsize: (optional; default: []) maximum size (cardinality) the cMCS
%      may have, use [] if you do not want to limit cMCS size
%
%   filename (optional): the calculated (c)MCSs can be stored in a file with the specified name.
%			If no filename is specified, results will not be stored.
%
%   reac_off (optional; default: []): a row vector specifying indices of reactions that
%      are considered to be inactive (i.e., their reaction rate is set to zero).
%
%   regulation (optional; default: []): is a structure that can be used to calculate regulatory cMCSs,
%      where reaction knockouts (cuts) can be combined with up- or down-regulaton
%      of certain reaction fluxes to solve the intervention problem. The following fields must be defined:
%
%	   regulation.reg_ind: vector with the indices of reactions whose flux can be regulated.
%
%          regulation.reg_down_up: a (2 x numel(regulation.reg_ind)) array indicating in the first row
%		which of the regulated reactions can be down-regulated (1) or not (0) and in the second
%		row which of the regulated reactions can be up_regulated (1) or not (0). A reaction may
%               be considered for both up- and downregulation. At least one '1' should appear for each
%               reaction (column) in this array.
%
%      In addition, one of the two fields must be specified:
%
%	   regulation.reg_bounds: is a cell array with numel(regulation.reg_ind) many entries storing
%               for each regulated reaction a vector with discrete flux levels to be considered when
%		up-/downregulating the reaction.
%
%          regulation.numregsteps: a number indicating how many different discrete flux levels are
%               to be considered. The boundaries are automatically determined as equidistant steps
%               in the feasible flux ranges of the regulated reactions.
%               Regulated reactions should therefore have bounded flux ranges for the desired space.
%               (regulation.numregsteps will be neglected if regulation.reg_bounds is given.)
%
%   For each flux level of each regulated reaction, an auxilliary reaction is introduced whose 'cut'
%   indicates that the regulated reaction is up/down-regulated with the respective level (see also the
%   returned arguments 'mcs' and 'reacNames').
%
%   time_limit: maximal time after which to stop calculation
%
%   EpGap: (optional; default: 0.98): relative MIP gap tolerance, use
%      smaller values to prolong the minimization phase
%
% Results:
% --------
%
%   mcs: the computed constrained minimal cut sets. Rows: cMCSs, columns: reactions.
%      A '1' in mcs(i,j) indicates the knockout of reaction j in mcs i. Note that the columns
%      correspond to the cnap.numr many reactions in cnap extended by the up/down-regulations
%      of the regulated fluxes (if regulation was considered). For the latter, a '1' means that
%      the respective regulation is active.
%
%   reacNames: a cell array with the (reaction) names of the columns in mcs. The first
%      cnap.numr entries correspond to the original reaction names given in cnap.reacID.
%      If regulated reactions where considered (see above), the names of these columns will
%      indicate the respective regulation, e.g. 'R2 >= 3.4' or 'R89 <= 48'.




mcs=[];
reg=[];
sys=[];
reacNames=[];
obj= [];
cs= [];

if(nargin<5)
    disp('This function needs at least 6 arguments');
    return;
end
if nargin < 6
    notknockable= [];
end
if nargin < 7
    maxMCSnumPre= 1;
end
if nargin < 8
    maxMCSnum= 0;
end
if nargin < 9
    maxMCSsize= [];
end
if(nargin<10)
    filename=[];
    dosave = 0;
elseif ~isempty(filename)
    dosave = 1;
else
    dosave = 0;
end
if(nargin<11)
    reac_off=[];
end


if(nargin<12 || isempty(regulation))
    reg_ind=[];
    reg_down_up=[];
    numregsteps=0;
    reg_bounds=[];
else
    if(~isfield(regulation,'reg_ind'))
        disp('Argument ''regulation'' requires mandatory field ''reg_ind'' .');
        return;
    else
        reg_ind=regulation.reg_ind;
    end
    
    if(isfield(regulation,'reg_down_up') && isfield(regulation,'reg_bounds') && ~isempty(regulation.reg_bounds))
        reg_down_up=regulation.reg_down_up;
        reg_bounds=regulation.reg_bounds;
        if(~iscell(reg_bounds) || numel(reg_bounds)~=numel(reg_ind))
            disp('Argument ''regulation.reg_bounds'' must be a cell array containing in the rows the vectors with the (regulated) flux bounds for each regulated variable.');
            return;
        elseif(size(reg_down_up,1)~=2 || size(reg_down_up,2)~=numel(reg_ind))
            disp('Argument ''reg_down_up'' must have two rows and numel(reg_ind) many columns.');
            disp('The first/second row indicates whether down/up-regulation is to be considered for the respective regulated reaction.');
            return;
        else
            for i=1:numel(reg_ind)
                reg_bounds{i}=sort(reg_bounds{i});
            end
            numregsteps=0;
        end
    elseif(isfield(regulation,'reg_down_up') && isfield(regulation,'numregsteps'))
        reg_bounds=[];
        reg_down_up=regulation.reg_down_up;
        numregsteps=regulation.numregsteps;
        if(size(reg_down_up,1)~=2 || size(reg_down_up,2)~=numel(reg_ind))
            disp('Argument ''reg_down_up'' must have two rows and numel(reg_ind) many columns.');
            disp('The first/second row indicates whether down/up-regulation is to be considered for the respective regulated reaction.');
            return;
        end
    else
        disp('Argument ''regulation'' requires field ''reg_bounds'' OR ''reg_down_up'' and ''numregsteps''.');
        return;
    end
end

if nargin < 13
    time_limit1=Inf;
end
if nargin < 14
    time_limit2=Inf;
end
if nargin < 15
    EpGap= 0.98;
elseif isempty(EpGap)
    EpGap= 0.98;
end
if nargin < 16
    inverseCutCount = [];
end

inhf=T;
ub=t;
des=D;
db=d;
no_cuts=notknockable;

st=initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
irr= cnap.reacMin >= 0;
n= size(st, 2);

zw=find(abs(cnap.reacMin)<=cnap.epsilon);
zw=zw(find(abs(cnap.reacMax(zw))<=cnap.epsilon));
reac_off=unique([reac_off,zw']);
cnap.reacMin(reac_off)=0;
cnap.reacMax(reac_off)=0;
if(~isempty(reac_off))
    disp('Reactions with specified zero flux:')
    disp(num2str(reac_off)) %P just indices (takes less space)
    % cnap.reacID(reac_off,:) %P full reaction names
end

cplex_inner= setup_cplex_inner_class_access();

javastderr= java.lang.System.err;
java.lang.System.setErr(java.io.PrintStream('cplex_stderr.log'));
rMin = cnap.reacMin;
rMin(rMin <= -1e3) = -Inf;
rMax = cnap.reacMax;
rMax(rMax >= 1e3) = Inf;
[sys,resok]= prepareMCSenum(st, reac_off, [no_cuts,reg_ind], inhf, ub, des, db,cnap.reacMin,cnap.reacMax,cnap.reacID,cnap.epsilon,inverseCutCount);
if(resok==0)
    return;
end
sys.des_rd_fvalb(sys.des_rd_fvalb==-Inf) = -1e3;
sys.des_rd_fvaub(sys.des_rd_fvaub==Inf) = 1e3;

%cnap.reacID(find(sys.blocked_fva),:)

java.lang.System.setErr(javastderr);
sys= add_subset_names(sys, cellstr(cnap.reacID));
sys.reg_ind=reg_ind;

[sub_reg_ind,dum]= find(sys.sub(:, reg_ind));
if length(unique(sub_reg_ind)) < length(reg_ind)
    zw=~any(sys.sub(:, reg_ind),1);
    disp(' ');
    disp('WARNING: the following regulated reactions are strictly coupled:');
    cnap.reacID(reg_ind(zw),:)
    disp('Select other reactions for regulaton. Exit.');
    disp(' ');
    return;
end


if(~isempty(reg_ind))
    disp('Feasible range of regulated reactions: full system LB/UB; reduced system LB/UB');
    disp([sys.des_fvalb(reg_ind), sys.des_fvaub(reg_ind), sys.des_rd_fvalb(sub_reg_ind), sys.des_rd_fvaub(sub_reg_ind)])
    
    interval_bounds= cell(length(sub_reg_ind), 1);
    if(isempty(reg_bounds))
        if(any(isinf(sys.des_rd_fvalb(sub_reg_ind))) || any(isinf(sys.des_rd_fvaub(sub_reg_ind))))
            disp('Some of the regulated reactions have unbounded flux ranges -> cannot determine suitable interval for flux regulation.')
            disp('Please ensure bounded fluxes for regulated reactions.');
            return;
        end
        
        num_intervals= (numregsteps+1)*ones(1, length(sub_reg_ind)); % must all be >= 2
        for i= 1:length(sub_reg_ind)
            interval_bounds{i}= linspace(sys.des_rd_fvalb(sub_reg_ind(i)), sys.des_rd_fvaub(sub_reg_ind(i)), num_intervals(i) + 1);
        end
    else
        for i= 1:length(sub_reg_ind)
            if(any(sys.des_rd_fvalb(sub_reg_ind(i))>reg_bounds{i}) || any(sys.des_rd_fvaub(sub_reg_ind(i))<reg_bounds{i}))
                disp(['Error: specified regulation levels for reaction ',sys.react_name{reg_ind(i)},' out of feasible range.']);
                fprintf('Regulation levels: %.2f \n', reg_bounds{i});
                fprintf('Feasible range: %.2f ... %.2f\n', sys.des_rd_fvalb(sub_reg_ind(i)),sys.des_rd_fvaub(sub_reg_ind(i)));
                return;
            end
            interval_bounds{i}= [sys.des_rd_fvalb(sub_reg_ind(i)), reg_bounds{i}, sys.des_rd_fvaub(sub_reg_ind(i))];
        end
    end
    
    for i= 1:length(sub_reg_ind)
        fprintf('%s\t', sys.react_name{reg_ind(i)});
        fprintf('%.2f ', full(sys.sub(sub_reg_ind(i), reg_ind(i)))*interval_bounds{i});
        fprintf('\n');
    end
    
    reg= add_flux_bound_interventions(sys, sub_reg_ind, interval_bounds,reg_down_up);
else
    reg=sys;
    reg.slack_groups=[];
end

allowed_cuts=reg.cuts;
allowed_cuts(any(reg.sub(:, no_cuts), 2))= false;

%options.cuts= reg.cuts;
%options.cuts(any(reg.sub(:, no_cuts), 2))= false;
%options.validate_mcs= true;

if (any(isinf(reg.des_rd_fvalb)) || any(isinf(reg.des_rd_fvaub)))
    disp('Problem: there are unbounded fluxes for desired system --> not allowed in integrated mode! ');
    disp('Please change flux bounds accordingly!');
    return;
end

%% integrated MILP
flux_lb= reg.des_rd_fvalb; %flux bounds for desired system!
flux_ub= reg.des_rd_fvaub; %flux bounds for desired system!

obj = ConstrainedMinimalCutSetsEnumerator2(reg.rd,reg.irrev_rd, [],reg.inh_rd,reg.ub,allowed_cuts,flux_lb,flux_ub,reg.des_rd,reg.db,reg.invCuts_rd);
obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 4096);
for i= 1:length(reg.slack_groups)
    obj.cpx.addLe(obj.cpx.sum([obj.z_vars(reg.slack_groups{i}), obj.z_vars(reg.slack_groups{i}+size(reg.rd, 2))]), 1);
end
disp(' ');
disp('Starting main procedure ...');
disp(' ');
tol= 1e-8;
obj= set_objective_function2(obj);
obj.cpx.setParam(cplex_inner.DoubleParam.EpGap, EpGap);
obj.cpx.setParam(cplex_inner.DoubleParam.EpInt, tol);
obj.cpx.setParam(cplex_inner.DoubleParam.EpRHS, tol);
obj.cpx.setParam(cplex_inner.DoubleParam.EpOpt, tol);
obj.cpx.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 4096);
if ~isempty(maxMCSsize)
    obj.evs_sz.setUB(maxMCSsize);
else
    maxMCSsize = length(reg.irrev_rd);
    obj.evs_sz.setUB(maxMCSsize);
end
if maxMCSnumPre == 0
    enum_method = 2;
elseif maxMCSnumPre == 1
    enum_method = 3;
elseif maxMCSnumPre > 1
    enum_method = 4;
else
    error('please set positive maxMCSnum');
end
if enum_method == 3 || enum_method == 4           %max_ev_size, num_evs_lim, time_limit, enum_method, validate, cplex_parset,inverseCutCount
    [obj, cs, status]= shortest_minimal_cut_sets2(obj,maxMCSsize,maxMCSnumPre,time_limit1,enum_method,0,[],reg.invCuts_rd);
    infeasible= status.equals(cplex_inner.Status.Infeasible);
    if ~isempty(cs) && ~infeasible && enum_method >= 3
        lplbrd= -Inf(length(reg.irrev_rd), 1);
        lplbrd(reg.irrev_rd~= 0)= 0;
        lpubrd= Inf(length(reg.irrev_rd), 1);
        [is_cs, dummy, is_des]= validate_mcs(reg.rd, reg.irrev_rd, reg.inh_rd, reg.ub, cs,...
            reg.des_rd, reg.db, lplbrd, lpubrd, 1);
        for i = find(is_cs & is_des) % if is_cs && is_des
            mcs2= make_minimal_cs(reg.rd, reg.irrev_rd, reg.inh_rd, reg.ub, cs(:,i), lplbrd, lpubrd);
            [is_cs2, is_min2, is_des2]= validate_mcs(reg.rd, reg.irrev_rd, reg.inh_rd, reg.ub,...
                mcs2, reg.des_rd, reg.db, lplbrd, lpubrd);
            if is_cs2 && is_min2 && is_des2
                cs(:,i)= mcs2;
            else
                disp('Could not validate cMCS, returning constrained cut set instead.');
            end
        end
    end
    if infeasible
        disp('No cMCS possible.');
    elseif isempty(cs)
        disp('No cMCS found.')
    end
end

%P Do calculation of cMCS starting from 1
%% copy-paste from CNAMCSEnumerator
if ~isempty(maxMCSnum) && (~isempty(cs) || enum_method == 2)
    reg.des_rd = reg.des2_rd;
    reg.db = reg.db2;
    q= size(reg.rd, 2);
    mcs = repmat(0,0,q);
    if enum_method ~= 2 %P in case of pre-calculation of big cut sets / find
        allowed_cuts = sum(cs,2) ~= 0;
        maxMCSsize = max(sum(cs));
    end
    obj= ConstrainedMinimalCutSetsEnumerator2(reg.rd,reg.irrev_rd, [],reg.inh_rd,reg.ub,allowed_cuts,[],[],[],[],reg.invCuts_rd);
    obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 4096);
    obj.cpx.setParam(cplex_inner.IntParam.Threads, 2);
    %Preparing tests for constrained cut sets
    if(~isempty(reg.des_rd))
        cgp=Cplex();
        cgp.Param.emphasis.numerical.Cur= 1;
        cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
        cgp.Param.simplex.tolerances.feasibility.Cur= cgp.Param.simplex.tolerances.feasibility.Min;
        cgp.DisplayFunc=[];
        cgp.Model.A= reg.rd;
        cgp.Model.ub= Inf(q, 1);
        cgp.Model.lb= -Inf(q, 1);
        cgp.Model.lb(reg.irrev_rd ~= 0)= 0;
        cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
        cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
        
        cgp.Model.A(end+1:end+size(reg.des_rd,1),:)=reg.des_rd;
        cgp.Model.lhs(end+1:end+size(reg.des_rd,1))=-Inf;
        cgp.Model.rhs(end+1:end+size(reg.des_rd,1))=reg.db;
        
        cgp.Model.obj=zeros(q, 1);
        cgp.Model.sense= 'maximize';
    end
    tot_comp_time = 0;
    for i=1:maxMCSsize
        if(dosave)
            fp=fopen(filename,'w');
            fclose(fp);
            obj.ev_size_lb = i;
            [obj, mcsn, status, comp_time]= shortest_minimal_cut_sets2(obj,i,maxMCSnum,time_limit2-tot_comp_time,2,[],[],reg.invCuts_rd); %cannot take maxMCSnum here since many will be filtered afterwards
            
            tot_comp_time = tot_comp_time+comp_time;
            
            if ~obj.err
                
                [mcsn,~,~] = unique(mcsn','rows','legacy'); % eliminate redundant mcs
                mcsn=mcsn';
                
                disp(' ');
                disp([num2str(size(mcsn,2)),' (compressed) MCSs of size ',num2str(i),' found.']);
                
                %% LP check to filter cMCS from MCS
                
                if(~isempty(mcsn))
                    if(~isempty(reg.des_rd))
                        disp(' ');
                        disp('Testing which MCSs are constrained MCSs ...');
                        ct= cputime;
                        [mcsn,~,~] = unique(mcsn','rows','legacy');
                        mcsn = mcsn';
                        [feasible, objmax]= lp_check_cmcs(mcsn, cgp);
                        disp(cputime - ct);
                        
                        disp(' ');
                        disp([num2str(sum(feasible)),' (compressed) constrained MCSs of size ',num2str(i),' found.']);
                        
                        newmcs=expand_mcs(mcsn(:, feasible), reg.sub)';
                        disp(' ');
                        disp([num2str(size(newmcs,1)),' (uncompressed) constrained MCSs of size ',num2str(i),' found.']);
                        disp(' ');
                        if(~isempty(newmcs))
                            mcs= [mcs;mcsn(:, feasible)'];
                        end
                    else
                        newmcs=expand_mcs(mcsn, reg.sub)';
                        disp([num2str(size(newmcs,1)),' (uncompressed) MCSs of size ',num2str(i),' found.']);
                        disp(' ');
                        if(~isempty(newmcs))
                            mcs= [mcs;mcsn'];
                        end
                    end
                    
                    fp=fopen(filename,'w');
                    fprintf(fp,['MCS up to size ',num2str(i),':\n']);
                    
                    for j=1:size(mcs,1)
                        zw=find(mcs(j,:));
                        for k=1:numel(zw)
                            fprintf(fp,[deblank(cnap.reacID(zw(k),:)),'	']);
                        end
                        fprintf(fp,'\n');
                    end
                    fclose(fp);
                end
                if size(mcs,1) >= maxMCSnum
                    break;
                end
            else
                mcs= [mcs;cs'];
                fprintf('There was an error during mcs enumeration\n');
                break;
            end
            disp(' ');
            disp([num2str(size(mcs,1)),' MCSs up to size ',num2str(i),' found so far.']);
            disp(' ');
            if(size(mcs,1)>=maxMCSnum)
                break;
            end
            % take 'false' away to get mails
            if false && strcmp(getpref('Internet','SMTP_Server'),'group.mpi-magdeburg.mpg.de') && i>=6
                [pth,fle,~]=fileparts(filename);
                sendmail('schneiderp@mpi-magdeburg.mpg.de',['Otto: ' fle ' cMCS size ' num2str(i) ' completed'],...
                    [[num2str(size(mcs,1)),' cMCSs up to size ',num2str(i),' found.'] 10 [num2str(floor(toc/60/60/24)) ' days ' num2str(floor(mod(toc/60/60,24))) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's']]);
                save([pth '/' fle '_cMCS' '.mat'], 'mcs', 'cnap');
            end
            
            disp(' ');
            if(isempty(D))
                disp(['Final result: ',num2str(size(mcs,1)),' MCSs found.']);
            else
                disp(['Final result: ',num2str(size(mcs,1)),' constrained  MCSs found.']);
            end
            
        else
            [obj, mcsn, status, comp_time]= shortest_minimal_cut_sets2(obj,maxMCSsize,maxMCSnum,time_limit2-tot_comp_time,2,[],[],reg.invCuts_rd);
            
            tot_comp_time = tot_comp_time+comp_time;
            
            if ~obj.err
                
                disp(' ');
                disp([num2str(size(mcsn,2)),' (compressed) MCSs of size ',num2str(i),' found.']);
                
                %% LP check to filter cMCS from MCS
                
                if(~isempty(reg.des_rd))
                    disp(' ');
                    disp('Testing which MCSs are constrained MCSs ...');
                    ct= cputime;
                    [mcsn,~,~] = unique(mcsn','rows','legacy');
                    mcsn = mcsn';
                    [feasible, objmax]= lp_check_cmcs(mcsn, cgp);
                    disp(cputime - ct);
                    
                    disp(' ');
                    disp([num2str(sum(feasible)),' (compressed) constrained MCSs of size ',num2str(i),' found.']);
                    
                    newmcs=expand_mcs(mcsn(:, feasible), reg.sub)';
                    disp(' ');
                    disp([num2str(size(newmcs,1)),' (uncompressed) constrained MCSs of size ',num2str(i),' found.']);
                    disp(' ');
                    if(~isempty(newmcs))
                        mcs= [mcs;mcsn(:, feasible)'];
                    end
                else
                    newmcs=expand_mcs(mcsn, reg.sub)';
                    disp([num2str(size(newmcs,1)),' (uncompressed) MCSs of size ',num2str(i),' found.']);
                    disp(' ');
                    if(~isempty(newmcs))
                        mcs= [mcs;mcsn'];
                    end
                end
            else
                mcs= [mcs;cs'];
                fprintf('There was an error during mcs enumeration \n');
            end
        end
    end
end
if ~isempty(mcs)
    cs = mcs';
end
%% expand mcs

if(~isempty(cs))
    mcs= expand_mcs(cs, reg.sub)';
    % mcs = abs(cs'*reg.sub); % old: mcs= expand_mcs(cs, reg.sub)';
    % sorting for size
    
    reacNames=cell(1);
    for i=1:n
        reacNames{i,1}=sys.react_name{i};
    end
    mcs = [mcs,sum(mcs')']; %P Put number of cuts in last column
    mcs = sortrows(mcs,size(mcs,2)); %P Sort for last column
    mcs = mcs(:,1:end-1); %P Remove last column again
    
    if(~isempty(reg_ind))
        todel=n+1:2:size(mcs,2)-1;
        toadd=n+2:2:size(mcs,2);
        mcs(:,todel)=[];
        for i=1:numel(toadd)
            reacNames{n+i,1}=reg.react_name{toadd(i)};
        end
    end
end

if(~isempty(mcs))
    if(~isempty(filename))
        fp=fopen(filename,'w');
        fprintf(fp,['MCS up to size ',num2str(max(sum(mcs~=0,2))),':\n']);
        for j=1:size(mcs,1)
            zw=find(mcs(j,:));
            for k=1:numel(zw)
                fprintf(fp,[reacNames{zw(k)},'     ']);
            end
            fprintf(fp,'\n');
        end
        fclose(fp);
    end
    
    disp(' ');
    disp(['Total number of (expanded) MCSs found: ',num2str(size(mcs,1)),'.']);
    disp(' ');
    
    if ~isempty(inverseCutCount)
        mcs(:,inverseCutCount) = mcs(:,inverseCutCount).*2-1;
    end
end

% Mark "negative cuts" / 1 if reaction is not added; -1 if it is added


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% auxiliary functions %%%%%

function reg= add_flux_bound_interventions(sys, sub_reg_ind, interval_bounds, lbubsel)
% adds the pseudometabolites and pseudoreactions to the reduced system
% no constraints are added for the first/last element of interval_bounds as
% these are assumed to be the lower/upper boundaries of the fluxe values
% lbubsel can be used to select whether <= and >= are to be added for each
% regulated reactions

if(isempty(sub_reg_ind))
    reg=sys;
    reg.slack_groups=cell(1,1);
    return;
end

if nargin < 4
    lbubsel= true(2, length(sub_reg_ind)); % upper row activates <=, lower row >=
end

num_intervals= cellfun('length', interval_bounds) - 1;
num_pseudo_met= sum(lbubsel, 1) * (num_intervals(:) - 1);
num_pseudo_reac= 2*num_pseudo_met;
reg= struct();
reg.rd= sys.rd;
reg.rd(end+num_pseudo_met, end+num_pseudo_reac)= 0; % preallocate space
reg.inh_rd= sys.inh_rd;
reg.inh_rd(:, end+num_pseudo_reac)= 0;
reg.des_rd= sys.des_rd;
reg.des_rd(:, end+num_pseudo_reac)= 0;
reg.st= sparse(sys.st);
reg.sub= sys.sub;
reg.sub(end+num_pseudo_reac, end+num_pseudo_reac)= 0;
reg.react_name= sys.react_name;
reg.cuts= sys.cuts;
reg.cuts(end+1:end+num_pseudo_reac)= false;
reg.irrev_rd= sys.irrev_rd;
reg.irrev_rd(end+1:end+num_pseudo_reac)= 1;
reg.ub= sys.ub;
reg.db= sys.db;
reg.sub_name= sys.sub_name;
%reg.pr_lb= zeros(num_pseudo_reac/4, 2); %alt Axel
%reg.pr_ub= zeros(num_pseudo_reac/4, 2); %alt Axel
reg.pr_lb=zeros(0,2);
reg.pr_ub=zeros(0,2);
reg.slack_groups= cell(1, length(sub_reg_ind)); % groups slack reactions
reg.des_rd_fvalb= sys.des_rd_fvalb;
reg.des_rd_fvalb(end+1:end+num_pseudo_reac)= 0;
reg.des_rd_fvaub= sys.des_rd_fvaub;
reg.des_rd_fvaub(end+1:end+num_pseudo_reac)= Inf;

pm_ind= size(sys.rd, 1) + 1;
pr_ind= size(sys.rd, 2) + 1;
inh_ind= size(sys.inh_rd, 1) + 1;
des_ind= size(sys.des_rd, 1) + 1;
pr_lb_ind= 1;
pr_ub_ind= 1;
st_pr_ind= size(sys.st, 2) + 1;
for i= 1:length(sub_reg_ind)
    sub_ind= sub_reg_ind(i);
    sub_irr= reg.irrev_rd(sub_ind);
    for j= 2:length(interval_bounds{i})-1
        if lbubsel(2, i)
            reg.rd(pm_ind, sub_ind)= 1; % add production of pseudometabolite to regulated reaction
            reg.rd(pm_ind, pr_ind)= -1; % pseudoreaction for >= that cannot be knocked out
            reg.inh_rd(inh_ind, pr_ind)= -1;
            reg.ub(inh_ind)= -interval_bounds{i}(j);
            reg.des_rd(des_ind, pr_ind)= -1;
            reg.db(des_ind,1)= -interval_bounds{i}(j);
            reg.pr_lb(pr_lb_ind, :)= [pr_ind, interval_bounds{i}(j)];
            reg.irrev_rd(pr_ind)= sub_irr; %Axel old
            reg.des_rd_fvalb(pr_ind)= interval_bounds{i}(j);	%boundaries also to be for desired system (for integrated MILP)
            reg.des_rd_fvaub(pr_ind)= reg.des_rd_fvaub(sub_ind);	%boundaries also to be for desired system (for integrated MILP)
            reg.sub(pr_ind, st_pr_ind)= 1;
            reg.sub_name{pr_ind}= ' ';
            pr_ind= pr_ind + 1;
            inh_ind= inh_ind + 1;
            des_ind= des_ind + 1;
            pr_lb_ind= pr_lb_ind + 1;
            st_pr_ind= st_pr_ind + 1;
            reg.rd(pm_ind, pr_ind)= 1; % slack reaction that can be knocked out
            reg.cuts(pr_ind)= true;
            %reg.irrev_rd(pr_ind)= sub_irr; %Axel old
            
            %boundaries also to be for desired system (for integrated MILP)
            if ~sub_irr
                reg.des_rd_fvaub(pr_ind)= abs(reg.des_rd_fvaub(sub_ind))+abs(reg.des_rd_fvalb(sub_ind));  %to be safe
            else
                reg.des_rd_fvaub(pr_ind)= abs(reg.des_rd_fvaub(sub_ind));
            end
            %    reg.des_rd_fvaub(pr_ind)= inf;
            %Axel old
            %if sub_irr
            %  reg.des_rd_fvalb(pr_ind)= 0;
            %  reg.des_rd_fvaub(pr_ind)= reg.des_rd_fvaub(sub_ind);
            %end
            reg.slack_groups{i}(end+1)= pr_ind;
            %reg.sub_name{pr_ind}= sprintf('%s(%d) >= %.2f', reg.sub_name{sub_ind}, sub_ind,interval_bounds{i}(j));
            reg.sub_name{pr_ind}= sprintf('%s >= %.2f', sys.react_name{sys.reg_ind(i)}, interval_bounds{i}(j));
            
            reg.sub(pr_ind, st_pr_ind)= 1;
            reg.react_name{st_pr_ind}= reg.sub_name{pr_ind};
            st_pr_ind= st_pr_ind + 1;
            pr_ind= pr_ind + 1;
            pm_ind= pm_ind + 1;
        end
        
        if lbubsel(1, i)
            reg.rd(pm_ind, sub_ind)= 1; % add production of pseudometabolite to regulated reaction
            reg.rd(pm_ind, pr_ind)= -1; % pseudoreaction for <= that cannot be knocked out
            reg.inh_rd(inh_ind, pr_ind)= 1;
            reg.ub(inh_ind)= interval_bounds{i}(j);
            reg.des_rd(des_ind, pr_ind)= 1;
            reg.db(des_ind,1)= interval_bounds{i}(j);
            reg.pr_ub(pr_ub_ind, :)= [pr_ind, interval_bounds{i}(j)];
            reg.irrev_rd(pr_ind)= sub_irr; %Axel old
            %boundaries also to be for desired system (for integrated MILP)
            reg.des_rd_fvalb(pr_ind)= reg.des_rd_fvalb(sub_ind);
            reg.des_rd_fvaub(pr_ind)= interval_bounds{i}(j);
            reg.sub(pr_ind, st_pr_ind)= 1;
            reg.sub_name{pr_ind}= ' ';
            pr_ind= pr_ind + 1;
            inh_ind= inh_ind + 1;
            des_ind= des_ind + 1;
            pr_ub_ind= pr_ub_ind + 1;
            st_pr_ind= st_pr_ind + 1;
            reg.rd(pm_ind, pr_ind)= -1; % slack reaction that can be knocked out
            reg.cuts(pr_ind)= true;
            %reg.irrev_rd(pr_ind)= sub_irr; %Axel old
            
            %boundaries also to be for desired system (for integrated MILP)
            if ~sub_irr
                reg.des_rd_fvaub(pr_ind)= abs(reg.des_rd_fvaub(sub_ind))+abs(reg.des_rd_fvalb(sub_ind));  %to be safe
            else
                reg.des_rd_fvaub(pr_ind)= abs(reg.des_rd_fvaub(sub_ind));
            end
            %    reg.des_rd_fvaub(pr_ind)= inf;
            %Axel old
            %if sub_irr
            %  reg.des_rd_fvalb(pr_ind)= 0;
            %  reg.des_rd_fvaub(pr_ind)= reg.des_rd_fvaub(sub_ind);
            %end
            reg.slack_groups{i}(end+1)= pr_ind;
            %reg.sub_name{pr_ind}= sprintf('%s(%d) <= %.2f', reg.sub_name{sub_ind}, sub_ind,interval_bounds{i}(j));
            reg.sub_name{pr_ind}= sprintf('%s <= %.2f', sys.react_name{sys.reg_ind(i)}, interval_bounds{i}(j));
            
            reg.sub(pr_ind, st_pr_ind)= 1;
            reg.react_name{st_pr_ind}= reg.sub_name{pr_ind};
            st_pr_ind= st_pr_ind + 1;
            pr_ind= pr_ind + 1;
            pm_ind= pm_ind + 1;
        end
    end
end

reg.reacMin_rd=sys.reacMin_rd;
reg.reacMin_rd(end+1:pr_ind-1)=reg.des_rd_fvalb(numel(sys.reacMin_rd)+1:pr_ind-1);
reg.reacMax_rd=sys.reacMax_rd;
reg.reacMax_rd(end+1:pr_ind-1)=reg.des_rd_fvaub(numel(sys.reacMax_rd)+1:pr_ind-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sys= add_subset_names(sys, react_name)

sys.sub_name= cell(size(sys.sub, 1), 1);
for i= 1:size(sys.sub, 1)
    sys.sub_name{i}= '{ ';
    for j= find(sys.sub(i, :))
        sys.sub_name{i}= [sys.sub_name{i}, sprintf('%g %s(%d) ', full(sys.sub(i, j)), react_name{j}, j)];
    end
    sys.sub_name{i}= [sys.sub_name{i}, '}'];
end
sys.react_name= react_name;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sys, resok] = prepareMCSenum(st, reac_off, keep_single, inh, ub, des, db,reacMin,reacMax,reacID,epsilon,inverseCutCount)
% disable: inh <= ub
% keep: des <= db

irr=(reacMin>=0)';
resok=1;

sys.st= sparse(st);
sys.irr= irr;
sys.reacMin=reacMin;
sys.reacMax=reacMax;
[m, n]= size(st);

if isempty(des)
    %des=zeros(1,n); % Only temporaily
    des=repmat(0,0,n);
    withdes=0;
    disp('No desired flux vectors defined.');
else
    withdes=1;
end


CplexEpOpt= 1e-9; % 1e-9 is also used as default in the Java functions
cplex_inner= setup_cplex_inner_class_access();
emphasizeAccuracy= cplex_inner.ParameterSet.constructor.newInstance([]);
emphasizeAccuracy.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
% emphasizeAccuracy.setParam(cplex_inner.BooleanParam.PreInd, false);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpOpt, CplexEpOpt);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpRHS, CplexEpOpt);
% emphasizeAccuracy.setParam(cplex_inner.IntParam.ScaInd, 1);
emphasizeAccuracy.setParam(cplex_inner.IntParam.AdvInd, 0);

%%%%Check feasibility of desired and undesired fluxes

cgp= Cplex();
cgp.Param.emphasis.numerical.Cur= 1;
cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
cgp.Model.A= st;
cgp.Model.ub= reacMax;
cgp.Model.lb= reacMin;
cgp.Model.lhs= zeros(m, 1);
cgp.Model.rhs= zeros(m, 1);
cgp.Model.obj=zeros(n, 1);
cgp.Model.sense= 'maximize';
cgp.DisplayFunc=[];
cgp.addRows(-inf(size(inh,1),1),inh,ub);

x= cgp.solve();
if x.status == 3
    disp('Target flux vectors infeasible!');
    fprintf('status %d %s \n', x.status,x.statusstring);
    disp('Exit due to infeasibility of target flux vectors - no cuts required to block target flux vectors!');
    resok=0;
    return;
else
    disp('Target flux vectors feasible!');
end

if(withdes)
    cgp= Cplex();
    cgp.Param.emphasis.numerical.Cur= 1;
    cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
    cgp.Model.A= st;
    cgp.Model.ub= reacMax;
    cgp.Model.lb= reacMin;
    cgp.Model.lhs= zeros(m, 1);
    cgp.Model.rhs= zeros(m, 1);
    cgp.Model.obj=zeros(n, 1);
    cgp.Model.sense= 'maximize';
    cgp.DisplayFunc=[];
    cgp.addRows(-inf(size(des,1),1),des,db);
    
    x= cgp.solve();
    if x.status == 3
        disp('Desired flux vectors infeasible!');
        fprintf('status %d %s \n', x.status,x.statusstring);
        disp('Exit due to infeasibility of desired flux vectors!');
        resok=0;
        return;
    else
        disp('Desired flux vectors feasible!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Running FVA to find blocked reactions.\n')
sys.blocked_fva= blocked_reactions_fva(st, reac_off,reacMin,reacMax);
%sys.blocked_fva= blocked_reactions_fva(st, reac_off,irr);

%sys.blocked_fva(:)=0;  %if errors occur in preprocessing

fprintf('%d reactions are explicitly disabled, FVA finds %d additional blocked reactions.\n',...
    length(reac_off), sum(sys.blocked_fva) - length(reac_off));

if withdes % FVA to find reaction bounds under desired conditions
    fprintf('Running FVA to determine reaction bounds under desired conditions.\n')
    lhs= [zeros(m, 1); -Inf(length(db), 1)];
    rhs= [zeros(m, 1); db];
    lpub= reacMax;
    lplb= reacMin;
    %lplb(irr ~= 0)= 0;
    lplb(sys.blocked_fva)= 0;
    lpub(sys.blocked_fva)= 0;
    [ind_i, ind_j, val]= find([st; des]);
    fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    sys.des_fvalb= fvaj(1);
    sys.des_fvaub= fvaj(2);
    if(any(isnan(sys.des_fvalb)) || any(isnan(sys.des_fvaub)))
        disp('Error in preprocessing (while determining flux bounds for desired conditions)!');
        resok=0;
        return;
    end
    sys.des_ess= fvaj(2) < -CplexEpOpt | fvaj(1) > CplexEpOpt; % essential reactions
    fprintf('FVA identified %d reactions as essential for desired conditions.\n', sum(sys.des_ess));
    
    if(any(isinf(sys.des_fvalb)) || any(isinf(sys.des_fvaub)))
        disp('There are unbounded fluxes for desired system!');
    end
else %Determine reaction rate ranges for "normal system"
    fprintf('Running FVA to determine reaction bounds.\n')
    lhs= [zeros(m, 1)];
    rhs= [zeros(m, 1)];
    lpub= reacMax;
    lplb= reacMin;
    lplb(sys.blocked_fva)= 0;
    lpub(sys.blocked_fva)= 0;
    [ind_i, ind_j, val]= find([st]);
    fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    sys.des_fvalb= fvaj(1);
    sys.des_fvaub= fvaj(2);
    if(any(isnan(sys.des_fvalb)) || any(isnan(sys.des_fvaub)))
        disp('Error in preprocessing (while determining flux bounds)!');
        resok=0;
        return;
    end
    sys.des_ess=[];
end

fprintf('Removing blocked reactions and combining reaction subsets (compression).\n')
% reduction with implicit removal of blocked reactions
sys.single_reac= any(inh, 1) | any(des, 1);
sys.single_reac(keep_single)= true;
sys.single_reac= find(sys.single_reac);

cnap1.stoichMat=st;
cnap1.reacMin=-ones(size(st,2),1);
cnap1.reacMin(find(irr))=0;
cnap1.epsilon=epsilon;
cnap1.reacID=reacID;
cnap1=CNAgenerateMFNetwork(cnap1);
[sys.rd, sys.irrev_rd, sys.sub,sys.red_met_int, xxx]= CNAcompressMFNetwork(cnap1,sys.single_reac,[],1,0,1,find(sys.blocked_fva),0);
sys.sub=sys.sub';  %Transposed version required here (compatible with metatool)

[m, n]= size(sys.rd);
sys.inh_rd= inh*sys.sub';
sys.des_rd= des*sys.sub';
sys.ub= ub;
sys.db= db;
sys.reacMin_rd=zeros(n,1);
sys.reacMax_rd=zeros(n,1);
for i=1:n
    subidx=find(sys.sub(i,:));
    sys.reacMin_rd(i)=max(reacMin(subidx));
    sys.reacMax_rd(i)=min(reacMax(subidx));
end
%%% New compression end

sys.cuts= true(1, n);
if withdes
    %fprintf('Determining essential reactions for desired conditions in reduced system.\n');
    lplb= sys.reacMin_rd;
    lpub= sys.reacMax_rd;
    lhs= [zeros(m, 1); -Inf(length(db), 1)];
    rhs= [zeros(m, 1); db];
    [ind_i, ind_j, val]= find([sys.rd; sys.des_rd]);
    % calculate essential reactions by checking whether each reaction is a
    % cut set that disables the desired conditions
    %  valj= CplexValidateMCS.validate(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, 0:n-1, 0:n-1, emphasizeFeasibility);
    %  sys.ess_rd= find(valj(1));
    
    fprintf('Running FVA to determine reaction bounds under desired conditions in reduced system.\n')
    fvajrd= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    sys.des_rd_fvalb= fvajrd(1);
    sys.des_rd_fvaub= fvajrd(2);
    zw= fvajrd(2) < -CplexEpOpt | fvajrd(1) > CplexEpOpt; % essential reactions
    sys.ess_rd=find(zw);
    
    if all(any(sys.sub(sys.ess_rd, :), 1) == sys.des_ess')
        fprintf('%d reaction subsets in the reduced system are essential for desired conditions.\n', length(sys.ess_rd))
        sys.cuts(sys.ess_rd)= false;
    else
        warning('prepareMCSenum:essentialReactionMismatch',...
            'Essential reactions identified in the full system do not correspond to essential reactions found in the reduced system.');
    end
    
    fprintf('Validation of reduction with FVA bounds for desired conditions:\n')
    dev= validate_reduction_fva_bounds(sys.sub, sys.des_fvalb, sys.des_fvaub,...
        sys.des_rd_fvalb, sys.des_rd_fvaub);
    fprintf('%g ', max(dev));
    fprintf('(these values should not be much larger than %g).\n', CplexEpOpt);
    
    %Uncertain whether this is correct: have to check reversibility under both desired AND undesired bahvior!
    %  ind= find(~sys.irrev_rd);
    %  sel= sys.des_rd_fvalb(ind) > -CplexEpOpt;
    %  if any(sel)
    %    fprintf('Restricting %d reactions to irreversible according to FVA in reduced system.\n', sum(sel));
    %    sys.irrev_rd(ind(sel))= true;
    %  end
    %  sel= sys.des_rd_fvaub(ind) < CplexEpOpt;
    %  if any(sel)
    %    fprintf('Restricting %d reactions to irreversible with changed direction according to FVA in reduced system.\n', sum(sel));
    %    ind= ind(sel);
    %    sys.irrev_rd(ind)= true;
    %    sys= change_reaction_direction_in_reduced_system2(sys, ind);
    %  end
    
else % determine instead bounds for whole system
    lplb= sys.reacMin_rd;
    lpub= sys.reacMax_rd;
    lhs= [zeros(m, 1)];
    rhs= [zeros(m, 1)];
    [ind_i, ind_j, val]= find([sys.rd]);
    
    fprintf('Running FVA to determine reaction bounds in reduced system.\n')
    fvajrd= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    sys.des_rd_fvalb= fvajrd(1);
    sys.des_rd_fvaub= fvajrd(2);
    if(any(isnan(sys.des_rd_fvalb)) || any(isnan(sys.des_rd_fvaub)))
        disp('Error in preprocessing (while determining flux bounds in reduced system)!');
        resok=0;
        return;
    end
end

%integrate flux bounds in inh system
inhidx=size(sys.inh_rd,1);
sys.des2_rd = sys.des_rd;
sys.db2     = sys.db;
desidx=size(sys.des2_rd,1);

for i=1:n
    if(sys.reacMin_rd(i)~=0 && ~isinf(sys.reacMin_rd(i)))
        inhidx=inhidx+1;
        sys.inh_rd(inhidx,i)=-1;
        sys.ub(inhidx)=-sys.reacMin_rd(i);
        if ~isempty(sys.des2_rd)
            desidx=desidx+1;
            sys.des2_rd(desidx,i)=-1;
            sys.db2(desidx)=-sys.reacMin_rd(i);
        end
    elseif(sys.reacMin_rd(i)==0 & sys.irrev_rd(i)~=1)
        disp('Irreversibility definition!');
    end
    if(~isinf(sys.reacMax_rd(i)))
        inhidx=inhidx+1;
        sys.inh_rd(inhidx,i)=1;
        sys.ub(inhidx)=sys.reacMax_rd(i);
        if ~isempty(sys.des2_rd)
            desidx=desidx+1;
            sys.des2_rd(desidx,i)=1;
            sys.db2(desidx)=sys.reacMax_rd(i);
        end
    end
end

sys.invCuts_rd = [];

if exist('inverseCutCount','var')
    if ~isempty(inverseCutCount)
        sys.invCuts_rd = find(sys.sub*iv(size(sys.sub,2),inverseCutCount));
    end
end