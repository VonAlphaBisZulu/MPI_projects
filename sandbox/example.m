% 05.07.2019
% Script Philipp/Sven
% Extending metabolic networks with genes and gene-rule-pseudometabolites
% Computing Minimal Gene Cut Sets
%
% To avoid errors check:
% 1. Are you using the updated version of CNAfluxVariability (min. 7th June 2018)?
% 2. Do you have CPLEX installed and working properly?
%

try parpool; catch, end % Start parallel Pool to speed up FVA

cnap = CNAloadNetwork({'iJOcore';1},1,1); % load network

%% Reaction indices
idx.growth = find(strcmp(cellstr(cnap.reacID),{'BIOMASS_Ec_iJO1366_core_53p95M'}));
idx.atpm   = find(strcmp(cellstr(cnap.reacID),{'ATPM'}));
idx.prod   = find(strcmp(cellstr(cnap.reacID),{'EX_etoh_e'}));
cnap = CNAsetGenericReactionData(cnap,idx.prod,'grRules','a1 and a2 or a3'); % add genes to make reaction gene-KI-able
idx.acEx   = find(strcmp(cellstr(cnap.reacID),{'EX_ac_e'}));
idx.glcEx  = find(strcmp(cellstr(cnap.reacID),{'EX_glc__D_e'}));
idx.glycEx = find(strcmp(cellstr(cnap.reacID),{'EX_glyc_e'}));
idx.o2     = find(strcmp(cellstr(cnap.reacID),{'EX_o2_e'}));

%% Set uptake bounds
cnap.reacMin(idx.acEx) = -30;
cnap.reacMax(idx.acEx) = 1000;
cnap.reacMin(idx.glycEx) = -15;
cnap.reacMax(idx.glycEx) = 0;
cnap.reacMin(idx.glcEx)  = -10;
cnap.reacMax(idx.glcEx)  = 0;

%% Target and Desired Regions 
minYield = 0.3; % carbon yield
Cp = 2; % number of carbon atoms in product
Cs = 6; % number of carbon atoms in glucose
Ca = 2; % ... in acetate
Cg = 3; % ... in glycerol
% T1: Product/Ac+Glc+Glyc when Acetate exchange is <= 0
T1 = full(sparse( [1          1              1              1            2           3          ], ...
                  [idx.prod   idx.glcEx      idx.acEx       idx.glycEx   idx.acEx    idx.atpm   ], ...
                  [Cp         minYield*Cs    minYield*Ca    minYield*Cg  1           -1         ],3,cnap.numr));
t1 =  [  0  ; 0 ; -3.15];
% T2: Only Product/Glc+Glyc yield when Acetate is >= 0
T2 = full(sparse( [1          1                             1            2           3          ], ...
                  [idx.prod   idx.glcEx                     idx.glycEx   idx.acEx    idx.atpm   ], ...
                  [Cp         minYield*Cs                   minYield*Cg  -1          -1         ],3,cnap.numr));
t2 =  [  0  ; 0 ; -3.15];
% Growth >= 0.05
D1 = full(sparse( [1          2 ], ...
                 [idx.growth idx.atpm], ...
                 [-1         -1],2,cnap.numr));
d1 = [-0.05 ; -3.15];
% ATPM >= 8
D2 = full(sparse( 1   , ...
                 [idx.atpm], ...
                 -1   ,1,cnap.numr));
d2 = -8;

% read gene rules
[cnap,enzymes,genes,gr_rules] = CNAgenerateGERassociation(cnap);

% KO-cost vectors
rkoCost                                  = nan(cnap.numr,1);
rkoCost(idx.o2)                          = 1; % Oxygen uptake is knockable
rkiCost                                  = nan(cnap.numr,1);
rkiCost([idx.acEx idx.glcEx idx.glycEx]) = 0.1;
gkoCost                                  = ones(length(genes),1);
gkoCost(strcmp(genes,'s0001'))           = nan;
gkiCost                                  = nan(length(genes),1);
gkiCost(strcmp(genes,{'a1'}))            = 0.09; % KI of EX etoh
gkiCost(strcmp(genes,{'a2'}))            = 0.13; % KI of EX etoh
gkiCost(strcmp(genes,{'a3'}))            = 0.17; % KI of EX etoh

max_num_interv  = 4.5;
time_limit      = inf;
max_solutions   = inf;

[mcs1, cnap1, cmp_mcs1, cmp_cnap1, mcs_idx1] = CNAgeneMultiImposedMCSFinder(cnap, {T1,T2} , {t1,t2} , {D1,D2} , {d1,d2} ,...
                                                rkoCost,rkiCost, ... koCost, kiCost
                                                max_num_interv,time_limit,max_solutions,...
                                                1,0, ... use_compression,use_bigM,
                                                1,gkoCost,gkiCost,[],1); % enum_method, gkoCost, gkiCost, gr_rules, verbose
                                            
% First compression step is not always helpful. Sometimes probem size increases after adding the
% genes. Therefore use_compression can be set to 0, 1 (only second compression), 2 (both compression
% steps)