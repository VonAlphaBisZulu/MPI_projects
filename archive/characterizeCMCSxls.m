function validMCS = characterizeCMCSxls( cnap, cnapCore, mcs, idx, filepath, gmcs,gidx,enzymes,ecoliGeneNames, G0, uncert )
% cnap -> CNA project
%
% cnapCore -> core network or empty
%
% mcs -> Matrix that contains all mcs to be validated
%
% idx ist struct with fields
% idx.prod -> Reaction index of product export reaction
% idx.subs -> Reaction index of substrate source
% idx.bm   -> Reaction index of biomass growth
% idx.atpm -> Reaction index of ATP Maintanance
% idx.o2   -> Reaction index of O2 uptake
%
% filepath -> File path with filename and extension to save results

% startcna if necessary
if ~exist('cnan','var')
    startcna(1)
end

%% prepare model and table
cnapBU      = cnap;
% mcs = unique(mcs, 'rows'); % Would theoretically be good
% for gene-cut-sets, multiple gene cut-sets could end up to the same
% reaction cut-set

cnapCoreBU = [];

% prepare reduced model if existing
if ~isempty(cnapCore)
   
    % map GS to core model
    rindex = find(contains(strcat('*',cellstr(cnapBU.reacID),'*'),strcat('*',cellstr(cnapCore.reacID),'*')));
    rmap = sparse(1:length(rindex),rindex,1);
    %rmap(end,cnap.numr) = 0;
    rmap(cnapCore.numr,cnap.numr) = 0;
    sindex = find(contains(strcat('*',cellstr(cnapBU.specID),'*'),strcat('*',cellstr(cnapCore.specID),'*')));
    smap = sparse(1:length(sindex),sindex,1);
    smap(end,cnap.nums) = 0;

    
    if ismember(idx.prod,rindex) % if export reaction is already in core
        idx.prodC = find(rmap(:,idx.prod));
        cnapCoreBU  = cnapCore;
    else % if not, add an export reaction for target metabolite
        metName = strrep(cnap.specID(find(cnap.stoichMat(:,idx.prod)),:),'_e','_c');
        metIdx = find(strcmp(cellstr(cnapCore.specID),cellstr(metName)));
        if ~isempty(metIdx)
            cnapCore = CNAaddReactionMFN(cnapCore,['EX_' metName],[metName ' = '],0,1000);
            idx.prodC = cnapCore.numr;
            cnapCoreBU  = cnapCore;
        else % if metabolite is not found delete model
            disp('metabolite not in Core model')
        end
    end
    
    idx.subsC = find(rmap(:,idx.subs));
    idx.bmC   = find(rmap(:,idx.bm));
    idx.atpmC = find(rmap(:,idx.atpm));
    idx.o2C   = find(rmap(:,idx.o2));
    
    if ~isempty(cnapCore)
        mcscomp = mcs*rmap';
        if size(mcscomp,2)<cnapCore.numr
            mcscomp(end,cnapCore.numr) = 0;
        end
        mcscomp = [zeros(1,cnapCore.numr);iv(cnapCore.numr,idx.o2C)';mcscomp];
    end
    
    % Tests in enhanced compressed Network
    
    [cnapCoreEBU reac spec] = buildCoreModel(cnap,sindex',rindex',cellstr(cnap.reacID(idx.prod,:)));
    idx.o2C   = find(rmap(:,idx.o2));
    mcscompE = mcs*sparse(1:cnapCoreEBU.numr,reac,1,cnapCoreEBU.numr,cnap.numr)';
    mcscompE = [zeros(1,cnapCoreEBU.numr);iv(cnapCoreEBU.numr,find(reac==idx.o2))';mcscompE];
end

%% Prepare MDF simulation parameters
if nargin >= 10
    if ~isempty(G0)
        Cmin = 1e-6*ones(cnap.nums,1);
        Cmax = 0.02*ones(cnap.nums,1);
        Cmax(findStrPos(cnap.specID,'co2_c')) = 1e-4;
        Cmin(findStrPos(cnap.specID,'glc__D_e')) = 1e-6;
        Cmax(findStrPos(cnap.specID,'glc__D_e')) = 0.055557;
        idx_atp = findStrPos(cnap.specID,'atp_c');
        idx_adp = findStrPos(cnap.specID,'adp_c');
        idx_amp = findStrPos(cnap.specID,'amp_c');
        idx_nad = findStrPos(cnap.specID,'nad_c');
        idx_nadh = findStrPos(cnap.specID,'nadh_c');
        idx_nadp = findStrPos(cnap.specID,'nadp_c');
        idx_nadph = findStrPos(cnap.specID,'nadph_c');
        fixed_ratios(1,1:3) = [idx_atp  idx_adp   10];
        fixed_ratios(2,1:3) = [idx_adp  idx_amp    1];
        fixed_ratios(3,1:3) = [idx_nad  idx_nadh  10];
        fixed_ratios(4,1:3) = [idx_nadp idx_nadph 10];
        RT = 8.31446*300;
        if size(G0,1) < cnap.numr
            G0(end:cnap.numr) = NaN;
        end

        if size(uncert,1) < cnap.numr
            uncert(end:cnap.numr) = NaN;
        end
    end
else
    G0 = [];
end

%% prepare table headers
varNames = {'# Cuts';'Cuts'; 'r P min'; 'Y P/S min';'rBM';'r P min at rBM max';'r ATP max';'r P min at r ATP max';...
            '%r P min of best';'%r BM';'%r ATP max';'factor r P at r ATP max';'diff r P at r ATP max to r P min';...
            '# Acc Spec'};
        
if  ~isempty(cnapCore)
    varNames = [varNames; '# Cuts core enhanced';'rBM (E) max';'rP(E) at BM max'; '# Cuts core original';'rComp BM max';'rP at BM max'];
end

if ~isempty(G0)
    varNames = [varNames; 'MDF pathway yield'];
end

% add empty set as first row to mcs to calculate yields in the vanilla
% model
mcs = [zeros(1,cnap.numr);iv(cnap.numr,idx.o2)';mcs];

%% loop through mcs
for i = 1:size(mcs,1)
    %% set up model
    cnap = cnapBU;
    numcuts(i) = sum(abs(mcs(i,[1:idx.o2-1,idx.o2+1:end]))); % number of cuts w/o o2
    cnap.reacMin(find(mcs(i,:))) = 0;
    cnap.reacMax(find(mcs(i,:))) = 0;
    if i>2
        rowNames = [rowNames; {['mcs ' num2str(i-2)]}];
        % get cuts
        cuts(i) = cellstr(strjoin(strcat(cellstr(cnap.reacID(find(mcs(i,:)),:)), {' ('}, cellstr(num2str(find(mcs(i,:))')), {') .'})'));
        disp([num2str(i-2) ' of ' num2str(size(mcs,1)-2)]);
    elseif i == 1
        cuts(i) = {'.'};
        disp('aerobic');
        rowNames = {'aerobic'};
    elseif i == 2
        cuts(i) = {'.'};
        disp('anaerobic');
        rowNames = [rowNames; {'anaerobic'}];
    else
        error('something went wrong')
    end
    
    %% perform calculations
    cnap.objFunc = iv(cnap.numr,idx.prod);
    [~,~,~,rPMin(i)] = CNAoptimizeFlux(cnap,[],[],2);
    [YPSMin(i),~,~,~] = CNAoptimizeYield(cnap,-iv(cnap.numr,idx.prod)',-iv(cnap.numr,idx.subs)',[],[],2);
    cnap.objFunc = -iv(cnap.numr,idx.bm);
    [fv,~,~,rBMMax(i)] = CNAoptimizeFlux(cnap,[],[],2);
    rProdAtrBMmax(i) = fv(idx.prod);
    cnap.objFunc = -iv(cnap.numr,idx.atpm);
    [fv,~,~,~] = CNAoptimizeFlux(cnap,[],[],2);
    ATPMMax(i) = fv(idx.atpm);
    ATPM = ATPMMax(i) * iv(cnap.numr,idx.atpm,nan);
    cnap.objFunc = iv(cnap.numr,idx.prod);
    [~,~,~,rPmin_ASMax(i)] = CNAoptimizeFlux(cnap,ATPM,[],2);
    
    prBMMax(i) = 100*rBMMax(i)/rBMMax(1);
    pfluxATPMax(i) = 100*ATPMMax(i)/ATPMMax(1);
    factor_rP_ASMax(i) = rPmin_ASMax(i)/rPMin(i);
    diff_rP_ASMax(i) = rPmin_ASMax(i)-rPMin(i);
    
    %% check how many metabolites are still accessible and which can disrupt
    %% the growth-coupling
    
    prod = regexp(cellstr(cnap.specID), '.*_c$', 'match'); 
    prod = find(~cellfun(@isempty,prod));

    numr = cnap.numr;
    nums = cnap.nums;

    % one run to get a preselection
    % make reaction boundaries infinit and all metabolite export reaction
    % bounds to 1, then maximize to see how many metabolites are accessible
    % with open bounds. The total number of accessible metabolites must be 
    % a subset of this.
    
    cpx = Cplex();
    cpx.DisplayFunc = [];
    cpx.Model.A  = cnap.stoichMat;
    reacsCytEx = full(-sparse(prod,1:length(prod),1));
    cpx.Model.A(min(prod):max(prod),numr+1:numr+length(prod)) = reacsCytEx;
    cpx.Model.lb = [-1./~(cnap.reacMin<0)+1; zeros(length(prod),1)];
    cpx.Model.ub = [ 1./~(cnap.reacMax>0)-1; ones(length(prod),1)];
    cpx.Model.lhs = zeros(1,nums);
    cpx.Model.rhs = zeros(1,nums);
    cpx.Model.obj = [zeros(numr,1); ones(length(prod),1)];
    cpx.Model.sense = 'maximize';
    cpx.solve();
    IdxProdCandidates = find(cpx.Solution.x(numr+1:numr+length(prod)));
    
    % second run to iterate over the remaining metabolites
    
    tic
    cpx = Cplex();
    cpx.DisplayFunc = [];
    reacsCytEx = full(-sparse(IdxProdCandidates,1:length(IdxProdCandidates),1));
    cpx.Model.A  = cnap.stoichMat;
    cpx.Model.A(min(IdxProdCandidates):max(IdxProdCandidates),numr+1:numr+length(IdxProdCandidates)) = reacsCytEx;
    cpx.Model.lb = [cnap.reacMin; zeros(length(IdxProdCandidates),1)];
    cpx.Model.ub = [cnap.reacMax; zeros(length(IdxProdCandidates),1)];
    cpx.Model.obj = [zeros(numr,1); zeros(length(IdxProdCandidates),1)];
    cpx.Model.lhs = zeros(1,nums);
    cpx.Model.rhs = zeros(1,nums);
    cpx.Model.sense = 'maximize';

    numMetAccessible(i)=0;

    for j = 1:length(IdxProdCandidates)
        if j>1
            cpx.Model.obj(numr+j-1) = 0; % last metabolite export to 0;
            cpx.Model.obj(numr+j-1) = 0; % last metabolite export to 0;
        end
        cpx.Model.ub(numr+j)    = 1;
        cpx.Model.obj(numr+j)   = 1;
        cpx.solve();
        if logical(cpx.Solution.x(numr+j)>cnap.epsilon)
            numMetAccessible(i) = numMetAccessible(i)+1;
        end
    end
    toc
    
    
    %% test in enhanced Core model with necessary reactions added (but not other pathways)
    
    cnapCoreE = cnapCoreEBU;
    numcutsCompE(i) = sum(abs(mcscompE(i,[1:find(reac==idx.o2)-1,find(reac==idx.o2)+1:end]))); % number of cuts w/o o2
    cnapCoreE.reacMin(find(mcscompE(i,:))) = 0;
    cnapCoreE.reacMax(find(mcscompE(i,:))) = 0;
    cnapCoreE.objFunc = -iv(cnapCoreEBU.numr, find(reac==idx.bm));
    [fv,CompFeasibleE(i),~,rCompEBMmax(i)] = CNAoptimizeFlux(cnapCoreE,[],[],2);
    rCompEPAtBMmax(i) = fv(find(reac==idx.prod));
       
    %% test in original core model, where only one sink reaction is added
    if  ~isempty(cnapCoreBU)
        cnapCore = cnapCoreBU;
        numcutsComp(i) = sum(abs(mcscomp(i,[1:idx.o2C-1,idx.o2C+1:end]))); % number of cuts w/o o2
        cnapCore.reacMin(find(mcscomp(i,:))) = 0;
        cnapCore.reacMax(find(mcscomp(i,:))) = 0;
        cnapCore.objFunc = -iv(cnapCore.numr,idx.bmC);
        [fv,CompFeasible(i),~,rCompBMmax(i)] = CNAoptimizeFlux(cnapCore,[],[],2);
        rCompPAtBMmax(i) = fv(idx.prodC);
    else
        numcutsComp(i) = nan;
        rCompBMmax(i) = nan;
        rCompPAtBMmax(i) = nan;
    end
    
    %% if MDF data are available: Test if Output reaction is in best
    if nargin >= 10
        if ~isempty(G0)
            [mdf, v, conc, dfs]= CNAcomputeOptMDFpathway(cnap, RT, [G0 uncert], Cmin, Cmax,[],[],fixed_ratios);
            vMDFyield(i) = -v(idx.prod)/v(idx.subs);
        end
    end
end

if ~isempty(i)
    % minyield in comparison to the best
    maxMinr = max(rPMin);
    
    prPMin = 100*rPMin/maxMinr;
    
    B = [num2cell(numcuts'),cuts',num2cell([roundDec(rPMin',2), roundDec(-YPSMin',2) , roundDec(-rBMMax',3),roundDec([rProdAtrBMmax',ATPMMax',rPmin_ASMax',...
        prPMin' , prBMMax', pfluxATPMax',factor_rP_ASMax',diff_rP_ASMax'],2)]),num2cell(roundDec(numMetAccessible'))];
    if ~isempty(cnapCore)
        B = [B,num2cell([numcutsCompE',roundDec(-rCompEBMmax',3),roundDec(rCompEPAtBMmax',2),numcutsComp',roundDec(-rCompBMmax',3),roundDec(rCompPAtBMmax',2)])];
    end
    if ~isempty(G0)
        B = [B,num2cell(roundDec(vMDFyield',2))];
    end
    B = [cellstr(cnap.reacID(idx.prod,:)), varNames';cellstr(rowNames),B];
    %% when mcs are gene-mcs, also note what gene-knockouts need to be done
    if nargin >= 9
        if ~isempty(gmcs)
            geneIDMap = [num2cell([enzymes(:).genes]); [enzymes(:).strGene]];
            gName = {};
            for i = 1:size(gidx,2)
                posgIdent = find([geneIDMap{1,:}]==gidx(i),1);
                gidentifier = geneIDMap{2,posgIdent};
                gName(i) = ecoliGeneNames(strcmp(ecoliGeneNames(:,1),gidentifier),2);
            end
            gKOs = {};
            for i = 1:size(gmcs)
                gKOs1   = strcat(gName(gmcs(i,:)),',');
                numgKOs(i) = length(gKOs1);
                gKOs(i) = {[gKOs1{:}]};
            end
            B = [B, ['numGeneCuts';{0};{0};num2cell(numgKOs)'], ['geneKnockouts';{''};{''};gKOs']];
        end
    end
    disp(B);
    %% save model to file
    cell2csv(filepath,B,char(9));
else
    warning(['no cut sets for ' char(cnap.reacID(idx.prod,:))]);
end
end

