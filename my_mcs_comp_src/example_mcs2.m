% This script reproduces the results from Table 2, scenario 4. 
%
% MCS computation for the strain design of a 2,3-butanediol production host
% using co-feeding and supporting high ATP maintenance rates
%
% We enumerate all Minimal Gene Cut Sets up to the size of 6 for the
% strongly growth coupled production of 2,3-butanediol from glucose and/or
% acetate and/or glycerol with E. coli. Compared to the third scenario, a 
% second desired region is introduced to make sure all found strain designs 
% support higher ATP demands (18 mM/gBDW/h). GPR rule compression and network 
% compression are enabled to speed up the MCS computation.
%
% % required files/models:
%   iMLcore.mat
%
% % important variables:
%   max_num_interv - defines the maximum number of possible gene cuts and
%                    substrate additions
%
% % process:
%   0) Start parallel pool to speed up FVAs
%   1) Setup model, add heterologous  reactions
%   2) Define Target and Desired regions for MCS computation
%   3) Run MCS computation
%   4) Output
%
% Correspondence: schneiderp@mpi-magdeburg.mpg.de
% -Jun 2020
%

%% 0) Starting CNA and Parallel pool (for faster FVA), defining compression setting
if ~exist('cnan','var')
    startcna(1)
end
max_solutions   = 10;
max_num_interv  = 20;

options.milp_solver     = 'matlab_cplex'; % alternative: 'java_cplex'; 
options.preproc_D_violations = 0;
options.mcs_search_mode = 1; % searches MCS one by one. If you change this to 2, make sure all knock-out or knock-in costs are integers or nan

% If parallel toolbox available, start parallel pool.
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate'))
    parpool(); % remove this line if MATLAB Parallel Toolbox is not available
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end
options.preproc_check_feas = false;

%% 1) Model setup
% load model
load('iMLcore.mat')

% 2,3-butanediol pathways
cnap = CNAaddSpeciesMFN(cnap,'actn_c',0,'Acetoin');
cnap = CNAaddSpeciesMFN(cnap,'diact_c',0,'Diacetyl');
cnap = CNAaddSpeciesMFN(cnap,'23bdo_c',0,'3-Hydroxybutan-2-one'); % 2,3 butanediol
% 1st pathway
cnap = CNAaddReactionMFN(cnap,'ACLDC','1 alac__S_c + 1 h_c = 1 co2_c + 1 actn_c' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;alsD;:;}//END_GENERIC_DATA',0,0,0,0);
% 2nd pathway
cnap = CNAaddReactionMFN(cnap,'ALOX','1 alac__S_c + 1 h_c + 0.5 o2_c = 1 co2_c + 1 h2o_c + 1 diact_c' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;butA;:;}//END_GENERIC_DATA',0,0,0,0);
cnap = CNAaddReactionMFN(cnap,'ACTD','1 h_c + 1 nadh_c + 1 diact_c = 1 nad_c + 1 actn_c' ,-1000,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);
% final step and export
cnap = CNAaddReactionMFN(cnap,'BTDD','1 h_c + 1 nadh_c + 1 actn_c = 1 nad_c + 1 23bdo_c' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;budC;:;}//END_GENERIC_DATA',0,0,0,0);
cnap = CNAaddReactionMFN(cnap,'EX_23bdo_e','1 23bdo_c =' ,0,1000,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);

% add alternative substrate supplies
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
cnap = CNAaddReactionMFN(cnap,'EX_ac_up_e','1 ac_e =' ,-30,0,0,nan,0,...
'//START_GENERIC_DATA{;:;deltaGR_0;#;num;#;NaN;:;uncertGR_0;#;num;#;NaN;:;geneProductAssociation;#;str;#;;:;}//END_GENERIC_DATA',0,0,0,0);
cnap.reacMax(ismember(cnap.reacID,{'EX_glyc_e'})) = 0;
cnap.reacMin(ismember(cnap.reacID,{'EX_glyc_e'})) = -20;
% cnap.reacMin(ismember(cnap.reacID,{'NADH16pp'})) = -1000;

%% 2) Define MCS setup
% reaction indices used in Target and Desired region
r23BDO_ex = find(strcmp(cellstr(cnap.reacID),'EX_23bdo_e'));
rGlc_up  = find(strcmp(cellstr(cnap.reacID),'EX_glc__D_e'));
rGlyc_up = find(strcmp(cellstr(cnap.reacID),'EX_glyc_e'));
rAc_up   = find(strcmp(cellstr(cnap.reacID),'EX_ac_up_e'));
rAc_ex   = find(strcmp(cellstr(cnap.reacID),'EX_ac_e'));
rATPM    = find(strcmp(cellstr(cnap.reacID),'ATPM'));
rBM      = find(~cellfun(@isempty,(regexp(cellstr(cnap.reacID),'BIOMASS_.*_core_.*'))));

% Target region - Yield is now referred to carbon uptake with 23bdo/glc as
% the reference. The strain design task is to enforce a yield of 30% compared
% to the maximum possible yield.
fixed_fluxes = nan(cnap.numr,1);
fixed_fluxes([rAc_up,rGlyc_up]) = 0;
Ymax_23bdo_per_glc = CNAoptimizeYield(cnap,full(sparse(1,r23BDO_ex,1,1,cnap.numr)),full(sparse(1,rGlc_up,-1,1,cnap.numr)),fixed_fluxes);
Ymax_c = Ymax_23bdo_per_glc/6*4; % carbon related yield
Y_thresh = Ymax_c * 0.3; % 30 % of the maximum carbon yield
disp(['Minimum carbon product yield threshold set to ' num2str(Y_thresh)]);
% T1: Under all circumstances the 2,3 BDO / glc+glyc yiled should exceed
%     the yield threshold
T1 = full(sparse( [1         1          1          ], ...
                  [r23BDO_ex rGlc_up    rGlyc_up   ], ...
                  [4         6*Y_thresh	3*Y_thresh ],1,cnap.numr));
t1 =  0;
% T2: If Acetate is not secreted, the 2,3 BDO / glc+glyc+ac yiled should exceed
%     the yield threshold
T2 = full(sparse( [1         1          1           1           2       ], ...
                  [r23BDO_ex rGlc_up    rGlyc_up    rAc_up      rAc_ex  ], ...
                  [4         6*Y_thresh	3*Y_thresh	2*Y_thresh  1       ],2,cnap.numr));
t2 =  [  0 ; 0 ];

% Desired regions: 
% D1: Biomass yield equivalent to r_BM > 0.05 h^-1 at glucose uptake rate of 
%     10 mmol/h/gBDW (equivalent biomass/carbon yield when grown on other substrates)
% D2: ATPM >= 18 mM/gBDW/h
Y_BM = 0.005; % Minimum Biomass Yield (referred to glucose / 6C)
D1 = full(sparse( [1         1          1          1       ], ...
                  [rBM       rGlc_up    rGlyc_up   rAc_up  ], ...
                  [-6        -6*Y_BM    -3*Y_BM    -2*Y_BM ],1,cnap.numr));
d1 = 0;
D2 = full(sparse( 1, rATPM, -1,1,cnap.numr));
d2 = -18;

% knockable reactions: O2 exchange as a potential knockout
rkoCost = nan(cnap.numr,1);
rkoCost(strcmp(cellstr(cnap.reacID),'EX_o2_e')) = 0.5; % (knockout cost 0.5)
% knockable genes
[~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
gkoCost = ones(length(genes),1);
gkoCost(ismember(genes,'spontanous')) = nan;% pseudo-gene that marks spontanous reactions is not knockable

% addable reactions: glucose, glycerol or acetate supply
rkiCost = nan(cnap.numr,1);
rkiCost([rGlc_up rGlyc_up]) = 1; % adding glucose or glycerol would cost 1
rkiCost(rAc_up) = 0.5;           % adding acetate is cheaper (0.5)
% addable genes: 2,3-butanediol pathway 1 or 2
gkiCost = nan(length(genes),1);
gkiCost(ismember(genes,{'alsD','butA','budC'})) = 1; % (addition cost per gene: 1)


%% 3) MCS Computation
[rmcs, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx] = CNAgeneMCSEnumerator2(cnap, ...
                                                    {T1 T2}, {t1 t2}, ...
                                                    {D1 D2}, {d1 d2}, ...
                                                    rkoCost,rkiCost, ...  reaction KO cost, reaction addition cost
                                                    max_solutions,max_num_interv, ...
                                                    gkoCost,gkiCost, ...  gene KO cost, gene addition cost
                                                    [],options,... gpr_rules, options
                                                    1); % verbose

% gmcs: The true minimal cut sets on the entire gene- and reaction network
% gcnap: The gene-reaction network
% rmcs: The cut sets mapped to the original network. Effects of gene-knockouts are traced back to
%       reaction knockouts

%% 4) Output
% Of the results, 5 gene MCS are selected and returned as text
gmcs_selection = gmcs(:,round(linspace(1,size(gmcs,2),5)));
ko_text = cell(1,size(gmcs_selection,2));
for i = 1:size(gmcs_selection,2)
    kis = find(~isnan(gmcs_selection(:,i)) & gmcs_selection(:,i) > 0);
    kos = find(~isnan(gmcs_selection(:,i)) & gmcs_selection(:,i) < 0);
    for j = 1:length(kis)
        ko_text{j,i} = ['+ ' strtrim(gcnap.reacID(kis(j),:))];
    end
    ko_text{length(kis)+1,i} = '';
    for j = length(kis)+2:length(kis)+length(kos)
        ko_text{j,i} = ['- ' strtrim(gcnap.reacID(kos(j-length(kis)-1),:))];
    end
end
disp(ko_text);