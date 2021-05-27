% This script reproduces the results from Table 3:
% It computes the smallest MCS for the weakly and directionally growth-coupled
% and substrate uptake production of 10 different products with E. coli.
% pGCP: production potential at max growth rate
% wGCP: production at max growth rate
% dGCP: prodution at all positive growth rates
% SUCP: prodcution in all flux states
% Computations are repeated 12 times and have a time limit of 2 hours.
%
% MCS computation 
%
% % Content:
%   0) Start parallel pool (if available) to speed up computation
%   1) Setup model (add pathways if necessary)
%   2) Define Target and Desired regions for MCS computation
%   3) Run MCS computation
%   4) Validate MCS
%   5) Plot results
%
% Correspondence: cellnetanalyzer@mpi-magdeburg.mpg.de
% -May 2021
%

%% User settings

% select coupling type:
% potential growth-coupling
% weak growth-coupling
% directional growth-coupling
% substrate uptake coupling
coupling = 'potential growth-coupling';

% select product
%  1: ethanol
%  2: lysine
%  3: glutamate
%  4: isobutanol
%  5: 1,4-BDO
%  6: 2,3-BDO
%  7: itaconic acid
%  8: isoprene
%  9: butane
% 10: methacrylic acid
% 11: resveratrol
% 12: bisabolene
% 13: tryptophane
% 14: octyl acetate
% 15: 4-hydroxycoumarin

product = 12;
idx.subsYieldFactor = -6;
idx.prodYieldFactor = 15;

% select number and time limit for computations
num_iter = 4; % 12 computations
options.milp_time_limit = 7200; % 2h

% choose whether a new instance of MATLAB should be started for each MCS computation.
% This measure ensures memory cleanup between the MCS runs. (usually not needed with Windows)
solve_in_new_process = 1;

options.mcs_search_mode = 1; % find any MCS
maxSolutions = 5;
maxCost = 25;
verbose = 1;
atpm = 1;
gene_mcs = 1;

% If runnning on a system with a SLURM workload manager:
% Use directory on internal memory to share data between the workers. 
% If job is running as a SLURM ARRAY, the compression switches (and also other
% parameters if indicated) are overwritten
if ~isempty(getenv('SLURM_ARRAY_TASK_ID')) % overwrite options if a SLURM array is used
    setups = generate_SLURM_codes();
    [product,coupling,maxCost,gene_mcs,atpm] = derive_options_from_SLURM_array(setups(str2double(getenv('SLURM_ARRAY_TASK_ID'))));
end

%% 0) Starting CNA and Parallel pool (for faster FVA), defining computation settings
if ~exist('cnan','var')
    startcna(1)
end

if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate')) && ~solve_in_new_process
    % start parpool and locate preferences-directory to tmp path
    prefdir = start_parallel_pool_on_SLURM_node();
% If running on local machine, start parallel pool and keep compression
% flags as defined above.
%
% On a local machine without a SLURM workload manager, but MATLAB parallel toolbox installed:
elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) && ~solve_in_new_process %#ok<DCRENAME>
    parpool();
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end

[product_rID,species,reactions] = load_pathway(product);
clear settings
settings.product = product_rID;
settings.coupling = coupling;
settings.gene_mcs = gene_mcs;
settings.atpm = atpm;
settings.maxCost = maxCost;
settings.searchMode = options.mcs_search_mode;
disp(settings);

if ~solve_in_new_process
    if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate'))
        % start parpool and locate preferences-directory to tmp path
        prefdir = start_parallel_pool_on_SLURM_node();
    % If running on local machine, start parallel pool and keep compression
    % flags as defined above.
    %
    % On a local machine without a SLURM workload manager, but MATLAB parallel toolbox installed:
    elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
           (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) %#ok<DCRENAME>
        parpool();
        wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
    end
end

%% 1) Model setup
% load model from file
load(which('iML1515.mat'));
load(which('iML1515geneNames.mat'));
cnap = CNAcobra2cna(iML1515,0);
cnap = block_non_standard_products(cnap);
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
load(which('core.mat'));

% % Uncomment this for computation in core network
% load(which('core.mat'));
% % Reduce model to core network
% cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
% cnap = CNAdeleteSpecies(cnap,find(~ismember(cellstr(cnap.specID),core_specs)),0);

% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% Load heterologous pathways if necessary
for spec = species
    cnap = CNAaddSpeciesMFN(cnap,spec.spec_id,0,spec.spec_name);
    cnap = CNAsetGenericSpeciesData(cnap,cnap.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    cnap = CNAaddReactionMFN(cnap, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    cnap = CNAsetGenericReactionData(cnap,cnap.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end

% Checking mass and and charge balances after pathway additions
check_mass_balance(cnap);

% No minimum ATP Maintenance
if ~atpm
	cnap.reacMin(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'ATPM')))) = 0;
    cnap.reacMax(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'EX_glc__D_e')))) = -10;
end

% replace bounds with inf
cnap.reacMin(cnap.reacMin==-1000) = -inf;
cnap.reacMax(cnap.reacMax== 1000) =  inf;

if gene_mcs
    % For gene MCS
    % lump several gene subunits
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'ATPS4rpp')),'geneProductAssociation','atpS*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH16pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH17pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'NADH18pp')),'geneProductAssociation','nuo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD2')),'geneProductAssociation','frd*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'FRD3')),'geneProductAssociation','frd*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'CYTBO3_4pp')),'geneProductAssociation','cyo*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'THD2pp')),'geneProductAssociation','pnt*');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'PDH')),'geneProductAssociation','ace* and lpd');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'AKGDH')),'geneProductAssociation','sucAB and lpd');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'SUCOAS')),'geneProductAssociation','sucCD');
    cnap = CNAsetGenericReactionData(cnap,find(strcmp(cellstr(cnap.reacID),'SUCDi')),'geneProductAssociation','sdh*'); % sdhA,B,C,D
    [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
    % All genes are knockable apart from pseudo-gene that marks spontanous reactions
    gkoCost = ones(length(genes),1);
    gkoCost(ismember(genes,'spontanous')) = nan;
    % Reactions are not knockable apart from O2 supply
    rkoCost = nan(cnap.numr,1);
    rkoCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
    [full_cnap, rmap] = CNAintegrateGPRrules(cnap);
else
    % For reaction MCS
    % knockables: All reactions with gene rules + O2 exchange as a potential knockout
    gpr = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    koCost = double(cellfun(@(x) ~isempty(x),gpr));
    notknock_gene = ~cellfun(@isempty,regexp(gpr,'(spontanous|phoE|ompF|ompN|ompC)','match'));
%     notknock_gene = ~cellfun(@isempty,regexp(gpr,'(spontanous|phoE|ompF|ompN|ompC|sufA|sufB|sufC|sufD|sufE|sufS)','match')); % maybe also iscA|iscS|iscU
    koCost(koCost==0 | notknock_gene) = nan;
    koCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 1;
    full_cnap = cnap;
    rmap = eye(cnap.numr);
end

%% 2) Define MCS setup
substrate_rID = 'EX_glc__D_e';
biomass_rID   = 'BIOMASS_Ec_iML1515_core_75p37M';

modules{1}.sense = 'desired';
modules{1}.type  = 'lin_constraints';
[modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

idx.subs = find(ismember(cellstr(cnap.reacID),substrate_rID));
idx.bm = find(ismember(cellstr(cnap.reacID),biomass_rID));
idx.prod = find(ismember(cellstr(cnap.reacID),product_rID));

%% compute thresholds:
% wGCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate should be ensured at maximal growth after 
%       the interventions.
% dGCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate devided by the 20% of the maximal growth rate
%       is the ratio between production and growth that should be attained
%       in all flux states.
% SUCP: At 20% of maximal growth rate, maximize production rate. 
%       30% of this rate, devided by the maximal substrate uptake rate
%       should be attained in all flux states.

disp('Compute production (yield) thresholds.');
% 1. compute max growth
cnap.objFunc(:) = 0;
cnap.objFunc(ismember(cnap.reacID,{biomass_rID})) = -1;
fv = CNAoptimizeFlux(cnap,[],[],2,-1);
r_bm_max20 = 0.2*fv(ismember(cnap.reacID,{biomass_rID}));
% 2. compute max production at 20% growth
cnap.objFunc(:) = 0;
cnap.objFunc(ismember(cnap.reacID,{product_rID})) = -1;
fv_fix = nan(cnap.numr,1);
fv_fix(ismember(cnap.reacID,{biomass_rID})) = r_bm_max20;
fv = CNAoptimizeFlux(cnap,fv_fix,[],2,-1);
r_p_20 = 0.2*fv(ismember(cnap.reacID,{product_rID}));
% 3. Thresholds for dGCP and SUCP
Y_PBM = r_p_20/r_bm_max20;
Y_PS  = r_p_20/-fv(ismember(cnap.reacID,{substrate_rID}));

comptime = nan(4*num_iter,1);
gmcs_tot = nan(full_cnap.numr,0);
rmcs_tot = nan(cnap.numr,0);
for coupling = {'potential growth-coupling' 'weak growth-coupling' 'directional growth-coupling' 'substrate uptake coupling'}

    clear('modules');
    % Desired fluxes for all coupling cases: growth >= 0.05/h
    modules{1}.sense = 'desired';
    modules{1}.type  = 'lin_constraints';
    [modules{1}.V(1,:),modules{1}.v(1,:)] = genV([{biomass_rID} {'>='} 0.05 ],cnap);

    disp(' ');
    disp('=============');
    disp(coupling{:});
    disp('=============');
    switch coupling{:}
        case 'potential growth-coupling'
            modules{1}.type  = 'bilev_w_constr';
            % Desired flux states with maximal growth and positive ethanol production
            [modules{1}.V(2,:),modules{1}.v(2,:)] = genV([{product_rID} {'>='} r_p_20 ],cnap);
            [modules{1}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0    ],cnap); % maximize biomass
        case 'weak growth-coupling'
            modules{2}.sense = 'target';
            modules{2}.type  = 'bilev_w_constr';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{product_rID} {'<='} r_p_20 ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID} {'>='} 0.05   ],cnap);
            [modules{2}.c(1,:),~]                 = genV([{biomass_rID} {'>='} 0      ],cnap); % maximize biomass
            modules_wg = modules;
        case 'directional growth-coupling'
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / ' biomass_rID]}   {'<='}  Y_PBM ],cnap);
            [modules{2}.V(2,:),modules{2}.v(2,:)] = genV([{biomass_rID}   {'>='}  0.01 ],cnap);
        case 'substrate uptake coupling'
            modules{2}.sense = 'target';
            modules{2}.type = 'lin_constraints';
            [modules{2}.V(1,:),modules{2}.v(1,:)] = genV([{[product_rID ' / -' substrate_rID]}    {'<='}  Y_PS],cnap);
        otherwise
            error(['Coupling mode ''' char(coupling_i) ''' not found.']);
    end

    for i = 1:num_iter
        if solve_in_new_process    
            displ(['Running MCS computation in separate thread: ' num2str(i) '/' num2str(num_iter)],verbose);
            if gene_mcs  %#ok<*UNRCH>
               [rmcs, gmcs, comptime(i), full_cnap] = MCS_enum_thread(gene_mcs,cnap,modules,...
                                                rkoCost,...
                                                maxSolutions,maxCost,...
                                                gkoCost,...
                                                options,verbose);
            else
                [rmcs, ~, comptime(i)] = MCS_enum_thread(gene_mcs,cnap,  modules,...
                                                koCost,...
                                                maxSolutions,maxCost,...
                                                [],...
                                                options,verbose);
            end
        else
            displ(['Running MCS computation: ' num2str(i) '/' num2str(num_iter)],verbose);
            tic;
            if gene_mcs  %#ok<*UNRCH>
                [rmcs, gmcs, full_cnap] = ...
                    CNAgeneMCSEnumerator3(cnap,modules,...
                        rkoCost,[],...
                        maxSolutions,maxCost,...
                        gkoCost,[],[],...
                        options,verbose);
            else
                 [gmcs, status] = CNAMCSEnumerator3(cnap,modules,...
                        koCost,[],...
                        maxSolutions,maxCost,...
                        options,verbose);
            end
            comptime(i) = toc;
        end
        if ~isempty(gmcs)
            rmcs_tot  = [rmcs_tot,rmcs];
            gmcs_tot  = [gmcs_tot,gmcs];
        end
    end
end

if gene_mcs
    reac_names = full_cnap.reacID;
    reac_names(full_cnap.rType ~= 'r',:) = reac_names(full_cnap.rType ~= 'r',[4:end '   ']);
else
    reac_names = cnap.reacID;
end

filename = [fileparts(which('compute_bisabolene.m')) filesep 'results' filesep product_rID];
    
%% Characterization and ranking of MCS
% Instead of the gene-MCS, their corresponding reaction-representations are analyzed.
% This is preferred, because the reaction-model is smaller and therefore analysis is 
% faster than in the GPR-extended model. Furthermore different gene-MCS can lead to 
% identical 'phenotypes' when translated to the reaction-model and by analyzing rMCS
% only a reduced, non-redundant set of reaction-MCS needs therefore to be considered.
if full(~all(all(isnan(gmcs_tot)))) % if mcs have been found
    disp('Characterizing mcs');
  % 5.1) Lump redundant MCS and create flux bounds for each mutant model
    rmcs_tot(isnan(rmcs_tot)) = -inf; % this step allows to apply 'unique' to remove duplicates
    [rmcs_tot,~,gmcs_rmcs_map] = unique(rmcs_tot','rows');
    rmcs_tot = rmcs_tot';
    rmcs_tot(rmcs_tot == -inf) = nan;
    MCS_mut_lb = repmat({cnap.reacMin},1,size(rmcs_tot,2));
    MCS_mut_ub = repmat({cnap.reacMax},1,size(rmcs_tot,2));
    for i = 1:size(rmcs_tot,2)
        MCS_mut_lb{i} = MCS_mut_lb{i}.*(rmcs_tot(:,i)==1 | rmcs_tot(:,i)==0);
        MCS_mut_lb{i}(isnan(MCS_mut_lb{i})) = 0;
        MCS_mut_ub{i} = MCS_mut_ub{i}.*(rmcs_tot(:,i)==1 | rmcs_tot(:,i)==0);
        MCS_mut_ub{i}(isnan(MCS_mut_ub{i})) = 0;
    end
  % 5.2) Set relevant indices [criterion 2-7] and prepare thermodynamic (MDF) parameters [criterion 9]
    % reaction indices
    T = readcell(which('thermo_iML1515_ph7.50.csv'));
    for i = 1:numel(T)
        if ismissing(T{i})
            T{i} = nan;
        end
    end
    pos = findStrPos(cnap.reacID,T(2:end,2));
    deltaGR_0_raw = cell2mat(T(2:end,3));
    deltaGR_0 = nan(cnap.numr,1);
    deltaGR_0(pos(~~pos)) = deltaGR_0_raw(~~pos);
    
    uncertGR_0_raw = cell2mat(T(2:end,4));
    uncertGR_0 = nan(cnap.numr,1);
    uncertGR_0(pos(~~pos)) = deltaGR_0_raw(~~pos);
    
    cnap = CNAsetGenericReactionData_with_array(cnap,'deltaGR_0',num2cell(deltaGR_0));
    cnap = CNAsetGenericReactionData_with_array(cnap,'uncertGR_0',num2cell(uncertGR_0));
    
    [idx,mdfParam] = relev_indc_and_mdf_Param(cnap,idx);

  % 5.3) Define core metabolism [criterion 8]
    % Add the new reactions also to the list of reactions that will be
    % considered "core" reactions in the final MCS characterization and ranking
    new_reacs = ismember(cellstr(cnap.reacID),{reactions.reac_id});
    reac_in_core_metabolism = ismember(cellstr(cnap.reacID),core_reacs);
    reac_in_core_metabolism(new_reacs) = 1;
    lbCore = cnap.reacMin;
    ubCore = cnap.reacMax;
    lbCore(~reac_in_core_metabolism) = 0;
    ubCore(~reac_in_core_metabolism) = 0;
  % 5.4) Costs for genetic interventions  [criterion 10]
    intvCost                  = full_cnap.mcs.kiCost;
    intvCost(isnan(intvCost)) = full_cnap.mcs.koCost(isnan(intvCost));
    intvCost(full_cnap.rType == 'g') = 1;
    gene_and_reac_names = cellstr(full_cnap.reacID);
    gene_and_reac_names(full_cnap.rType == 'g') = cellstr(full_cnap.reacID((full_cnap.rType == 'g'),4:end)); % to avoid the 'GR-' prefix
    gene_and_reac_names(full_cnap.rType == 'g') = strrep(strrep(gene_and_reac_names(full_cnap.rType == 'g'),'(',''),')',''); % remove brackets
  % 5.5) Start characterization and ranking
    [MCS_rankingStruct, MCS_rankingTable]...
        = CNAcharacterizeGeneMCS( cnap , MCS_mut_lb, MCS_mut_ub, 1:size(MCS_mut_lb,2),... model, mutants LB,UB, incices ranked mcs
        idx, idx.cytMet, modules, mdfParam, ... relevant indices, Desired and Target regions
        lbCore, ubCore, gmcs_tot, intvCost, gene_and_reac_names, gmcs_rmcs_map, ...
        0:11, ones(1,11),2); % assessed criteria and weighting factors
    % save ranking and textual gmcs as tab-separated-values
    cell2csv([filename '.tsv'],MCS_rankingTable,char(9));
    text_gmcs = cell(size(gmcs_tot,2),1);
    for i = 1:size(gmcs_tot,2)
        mcs_tx = strtrim(strsplit(mcs2text(reac_names,gmcs_tot(:,i)),','));
        text_gmcs(i,1:numel(mcs_tx)) = mcs_tx;
    end
    cell2csv([filename '-gmcs.tsv'],text_gmcs,char(9));
    save([filename '.mat'],'MCS_rankingStruct','MCS_rankingTable','rmcs','gmcs_rmcs_map');
    if ~isempty(getenv('SLURM_JOB_ID'))
        system(['~/bin/sshpass -f ~/.kdas scp ' [filename '-gmcs.tsv'] ' schneiderp@linssh.mpi-magdeburg.mpg.de:/data/bio/teams/modeling/SchneiderP/Results/2020_mechthild']);
        system(['~/bin/sshpass -f ~/.kdas scp ' [filename '.mat'] ' schneiderp@linssh.mpi-magdeburg.mpg.de:/data/bio/teams/modeling/SchneiderP/Results/2020_mechthild']);
        system(['~/bin/sshpass -f ~/.kdas scp ' [filename '.tsv'] ' schneiderp@linssh.mpi-magdeburg.mpg.de:/data/bio/teams/modeling/SchneiderP/Results/2020_mechthild']);
    end
end
disp('Finished.');


%% Supplementary function. Generate Target or desired region from text.

function [V,v] = genV(constraints,cnap)
% Generate Vectors V and v so that V*r <= v
% input: some constraint seperated 3 three cells. e.g.:
%           r_1 + r_4 / r_3 - r_2    |    >=    |    a
    V = [];
    v = [];
    rMin = nan(cnap.numr,1);
    rMax = nan(cnap.numr,1);
    reacID = cellstr(cnap.reacID);

    for j = 1:size(constraints,1)
        % get right hand side
        a = constraints{j,3};
        if ischar(a)
            a = str2double(a);
            if isnan(a)
                error('error in ''t'' of target region definition');
            end
        end
        % get direction of inequality
        switch constraints{j,2}
            case '<='
                eqop = 1;
            case '>='
                eqop = -1;
            otherwise
                error('please define inequality either as ''<='' or ''>=''');
        end
        % split into numerator an divisor and get coefficients
        numDiv = strtrim(split(constraints{j,1},'/'));
        % get variables and coefficients for vector
        [num, cnum] = findReacAndCoeff(numDiv(1),reacID);
            % fractional constraint ispreprocessed
        if length(numDiv) == 2
            [div, cdiv] = findReacAndCoeff(numDiv(2),reacID);
            [rMin(div), rMax(div)] = CNAfluxVariability(cnap,[],[],-1,div,[],[],0);
            % check if all reactions take identical signs (+ or -)
            % this is needed to do the equation rearrangement. If the signs
            % are ambigous, a lot of case differentiations would be needed,
            % what is not done here.
            if any(rMin(div).*rMax(div) < 0)
                error(['reactions that are part of the divisor inside the '... 
                        'target constraints must not span a positive AND negative range']);
            end
            rDir = sign(rMin(div)+rMax(div));
            if any(rDir.*cdiv' > 0) &&  any(rDir.*cdiv' < 0)
                error(['reactions that are part of the divisor inside the '...
                       'target constraints must all have the same direction '...
                       '(positive or negative). Please check if coefficients and '...
                       'reaction ranges lead to all positive or all negative variables.']);
            end
            if sign(sum(rDir.*cdiv')) == -1 % if the divisor is all negative
                % change direction of inequality
                eqop = -eqop;
            end
            cdiv = -a*cdiv; % transformation to following form:
            a = 0;
            % (cnum num) - a*cdiv1*div1 - a*cdiv2*div2 <= 0
            num  = [num,   div];
            cnum = [cnum, cdiv];
        end
        % constraint is generated
        if eqop == -1  % if inequality is >= the whole term is multiplied with -1
            cnum = -cnum;
            a    = -a;
        end
        v(j,1)   = a;
        V(j,:) = full(sparse(num,1,cnum,length(reacID),1));
    end
end

function [ridx,coeff] = findReacAndCoeff(eq,reacID)
    coeff = [];
    r = cell.empty(1,0);
    ridx = [];
    for strPart = strsplit(char(eq))
        str = char(regexprep(strPart,'^(\s|-|\.|\()*|(\s|-|\.|\))*$','')); % remove leading and tailing special characters
        if ~isempty(str)
            v = regexp(reacID, ['^' str '$'], 'match');
            if any(~cellfun(@isempty,v))
                r(end+1) = {str};
                ridx = [ridx find(~cellfun(@isempty,v))'];
            end
        end
    end
    for k = 1:length(r)
        c = regexp(eq, ['(\s|\d|-|\.)*?(?=' r{k} '(\s|$))'], 'match');
        c = regexprep(char(c{:}),'\s','');
        switch c
            case ''
                coeff(k) = 1;
            case '+'
                coeff(k) = 1;
            case '-'
                coeff(k) = -1;
            otherwise
                coeff(k) = str2double(c);
        end
    end
end

% read certain cells
function C = readCells(table,keyword,cols) % returns all cells below the given keyword until the first empty line
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),keyword)));
    for i = 1:length(sheet)
        lastrow = row(i)+find(strcmp(strtrim(table(row(i):end,col(i),sheet)),''),1,'first')-2;
        if isempty(lastrow)
            lastrow = size(table,1);
        end
        C{i} = table((row(i)+1):lastrow,col(i):(col(i)+cols-1),sheet(i));
    end
    if ~exist('C','var')
        error(['ERROR: Keyword ''' keyword ''' was not found in any sheet.'])
    elseif length(C) == 1
        C = C{:};
    else
        disp(['WARNING: Keyword ''' keyword ''' was found multiple times in the sheets.']);
    end
end

function tab = readXls(filename)
    sheets = sheetnames(filename);
    for i=1:length(sheets)
        workSheet = readcell(filename,'sheet',i);
        if ~isempty(workSheet)
            tab(1:size(workSheet,1),1:size(workSheet,2),i) = workSheet;
        end
    end
    tab = strTable(tab);
    
    function rawmat = strTable(rawmat)
        [rows, cols, shts] = size(rawmat);
        for row = 1:rows
            for col = 1:cols
                for sh = 1:shts
                    if isnumeric(rawmat{row,col,sh})
                        if isempty(rawmat{row,col,sh}) || isnan(rawmat{row,col,sh})
                            rawmat(row,col,sh) = {''};
                        else
                            rawmat{row,col,sh} = num2str(rawmat{row,col,sh});
                        end
                    elseif isa(rawmat{row,col,sh},'missing')
                        rawmat(row,col,sh) = {''};
                    end
                end
            end
        end
    end
end

function [valid_T, valid_D] = verify_mcs(cnap,mcs,T,t,c,D,d)
    if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
           (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && ~isempty(gcp('nocreate')) %#ok<DCRENAME>
        numworkers = getfield(gcp('nocreate'),'NumWorkers');
    else
        numworkers = 0;
    end
    mcs = mcs<0;
    valid_T = nan(size(mcs,2),1);
    valid_D = nan(size(mcs,2),1);
    parfor (i = 1:size(mcs,2),numworkers)
        cnap_valid = cnap;
        cnap_valid.reacMin(mcs(:,i)) = 0;
        cnap_valid.reacMax(mcs(:,i)) = 0;
        cnap_valid_opt = cnap_valid;
        if ~isempty(c)
            cnap_valid.objFunc = c';
            fv = CNAoptimizeFlux(cnap_valid,[], [], 2, -1);
            cnap_valid_opt.reacMin(logical(c)) = fv(logical(c));
            cnap_valid_opt.reacMax(logical(c)) = fv(logical(c));
        end
        if isempty(T)
            valid_T(i) = 0;
            valid_D(i) = testRegionFeas(cnap_valid_opt,D,d,2);
        else
            valid_T(i) = testRegionFeas(cnap_valid_opt,T,t,2);
            valid_D(i) = testRegionFeas(cnap_valid,D,d,2);
        end
    end
end

%% Supplementary function. Only used in systems with SLURM workload management.
% digits 1-4  : product
% digits 5-7 : weak coupling, directional coupling, substrate coupling
% digits 8-16: maxCost
% digits 17: geneMCS
% digits 18: atpm
function [productID,coupling,maxCost,geneMCS,atpm] = derive_options_from_SLURM_array(numcode)
    settings = dec2bin(numcode,18);
    productID = bin2dec(settings(1:4));
    switch settings([5 6 7])
        case '000'
            coupling = 'auto';
        case '001'
            coupling = 'potential growth-coupling';
        case '010'
            coupling = 'weak growth-coupling';
        case '011'
            coupling = 'directional growth-coupling';
        case '100'
            coupling = 'ATP coupling';
        case '101'
            coupling = 'substrate uptake coupling';
    end
    maxCost = bin2dec(settings(8:16));
    geneMCS = bin2dec(settings(17));
    atpm = bin2dec(settings(18));
end

function setups = generate_SLURM_codes()
product = 1:15;
coupling = {'potential','weak','directional','substrate'};
maxCost = 60;
table = {}; % table only required for excel-file
for prod = product
    table(1+(prod-1)*5,1:4) = num2cell(4*(prod-1)+(1:4));
    a(1) = getSlurmArrayCode(product(prod),coupling{1},maxCost,1,1);
    a(2) = getSlurmArrayCode(product(prod),coupling{2},maxCost,1,1);
    a(3) = getSlurmArrayCode(product(prod),coupling{3},maxCost,1,1);
    a(4) = getSlurmArrayCode(product(prod),coupling{4},maxCost,1,1);
    table(2+(prod-1)*5,1:4) = num2cell(a);
    setups((prod-1)*4+(1:4)) = a;
end
end

function [rmcs, gmcs, comptime, full_cnap] = MCS_enum_thread(gene_mcs,cnap,modules,...
                                            rkoCost,...
                                            maxSolutions,maxCost,...
                                            gkoCost,...
                                            options,verbose)
        % start MCS Enumeration in separate MATLAB instance
        if ~isempty(getenv('SLURM_TEMP'))
            tempdir = getenv('SLURM_TEMP');
        else
            tempdir = eval('tempdir');
        end
        rng('shuffle');
        wdir = [tempdir num2str(randi([0,1e6-1])) filesep];
        mkdir(wdir);
        filename = [wdir 'ws_' num2str(randi([0,1e6-1])) '.mat'];
        save(filename,'gene_mcs','cnap','modules','rkoCost','maxSolutions','maxCost','gkoCost','options','verbose');
        wd = pwd;
        cd(evalin('base','cnan.cnapath'));
        system(['matlab -nosplash -nodesktop -r "addpath(''' genpath(fileparts(mfilename('fullpath'))) ''');' ... % add project path
               'addpath(''' evalin('base','cnan.cnapath') ''');startcna(1);' ... % add CNA path and start CNA
               'compute_mcs_ext(''' wdir ''',''' filename ''');exit()"']); % compute
        cd(wd);
        load(filename,'gmcs','rmcs','comptime','status','full_cnap');
        if status ~= 0
            comptime = nan;
        end
        rmdir(wdir,'s');
end

function indices = findStrPos( str , pattern , opts )
% str       the space that is seached
% pattern   the search keyword or pattern
% opts      options: 0 - normal search ||| 'regex' - regex search
    if nargin == 2
        opts = 0;
    end
    % convert str type to cell
    switch class(str)
        case 'string'
            str = cellstr(strtrim(str));
        case 'char'
            str = cellstr(str);
        case 'cell'

        otherwise
            error('input 1 of findStrPos doesn''t correct type');
    end
    % convert pattern type to cell
    switch class(pattern)
        case 'string'
            pattern = char(pattern);
        case 'char'
            pattern = cellstr(pattern);
        case 'cell'

        otherwise
            error('input 2 of findStrPos doesn''t correct type');
    end
    [rows,cols] = size(pattern);

    switch opts
        case 0
            for i = 1:(rows*cols)
                ind                      = find(strcmp(str, pattern(i)));
                indices(1:length(ind),i) = ind;
            end
        case 'regex'
            if opts
                for i = 1:(rows*cols)
                    match = regexp(str, pattern(i), 'match');
                    ind = find(~cellfun(@isempty,match));
                    indices(1:length(ind),i) = ind;
                end
            end
        otherwise
            error('define correct option. 0 for basic search or ''regex''');
    end
    if any( size(pattern) == 0) || (any(size(indices)==0))
        indices = [];
    end
end

function ko_ki_text = mcs2text(reac_names,mcs)
    reac_names = strrep(cellstr(reac_names),'_','\_');
    kos = find(mcs<0);
    kis = find(mcs>0);
    if ~isempty(kis) 
        kis = [strjoin(strcat('+',reac_names(kis)),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos) 
        kos = strjoin(reac_names(kos,:),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end

function [idx,mdfParam] = relev_indc_and_mdf_Param(cnap,idx)
% function is used to find reaction and species indices that are used for
% the characterization and ranking of MCS
    % relevant reaction indices
    idx.o2      = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*EX_o2_e.*','match')));
    idx.atpm    = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*ATPM','match')));
    % relevant species indices
    idx.cytMet = find(~cellfun(@isempty ,regexp(cellstr(cnap.specID),'_c$','match')));
    % other important species
    idx.pi    = find(strcmp(cellstr(cnap.specID), 'pi_c')); 
% (adding idx.pi activates the output of coupling mechanism analysis (not yet functional))
    idx.h     = find(strcmp(cellstr(cnap.specID), 'h_c'));
    idx.h2o   = find(strcmp(cellstr(cnap.specID), 'h2o_c'));
    idx.atp   = find(strcmp(cellstr(cnap.specID), 'atp_c'));
    idx.adp   = find(strcmp(cellstr(cnap.specID), 'adp_c'));
    idx.amp   = find(strcmp(cellstr(cnap.specID), 'amp_c'));
    idx.nad   = find(strcmp(cellstr(cnap.specID), 'nad_c'));
    idx.nadh  = find(strcmp(cellstr(cnap.specID), 'nadh_c'));
    idx.nadp  = find(strcmp(cellstr(cnap.specID), 'nadp_c'));
    idx.nadph = find(strcmp(cellstr(cnap.specID), 'nadph_c'));
    idx.co2_e = find(strcmp(cellstr(cnap.specID), 'co2_c'));
    idx.glc_e = find(strcmp(cellstr(cnap.specID), 'glc__D_e'));
    % MDF setup (thermodynamic benchmark)
    mdfParam.Cmin    = 1e-6*ones(cnap.nums,1);
    mdfParam.Cmin(idx.glc_e) = 1e-6;
    mdfParam.Cmax = 0.02*ones(cnap.nums,1);
    mdfParam.Cmax(idx.co2_e) = 1e-4;
    mdfParam.Cmax(idx.glc_e) = 0.055557;
    mdfParam.fixed_ratios(1,1:3) = [idx.atp   idx.adp   10];
    mdfParam.fixed_ratios(2,1:3) = [idx.adp   idx.amp    1];
    mdfParam.fixed_ratios(3,1:3) = [idx.nad   idx.nadh  10];
    mdfParam.fixed_ratios(4,1:3) = [idx.nadph idx.nadp  10];
    mdfParam.RT = 8.31446*300/1000; % Computation of MDF in kJ
    mdfParam.bottlenecks = 0; % change to 1 to compute thermodynamic bottlenecks
    mdfParam.G0 = cell2mat(CNAgetGenericReactionData_as_array(cnap,'deltaGR_0'));
    mdfParam.uncert = cell2mat(CNAgetGenericReactionData_as_array(cnap,'uncertGR_0'));
end
