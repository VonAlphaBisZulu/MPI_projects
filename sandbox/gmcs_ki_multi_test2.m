% 05.07.2019
% Philipp
% Testing the use of multiple substrates and multiple target and desired regions in iJOcore with gene MCS:
%
%% Preparation
cnap = CNAloadNetwork({'iJOcore';1},1,1);

cnap = CNAaddReactionMFN(cnap,'EX_ac_e_2'    ,'1 ac_e = '           ,-30,0);
cnap = CNAaddReactionMFN(cnap,'EX_glyc_e_2'  ,'1 glyc_e = '         ,-20,0);
cnap = CNAaddReactionMFN(cnap,'EX_succ_e_2'  ,'1 succ_e = '         ,-15,0);
cnap = CNAaddReactionMFN(cnap,'EX_glc__D_e_2','1 glc__D_e = '       ,-10,0);

idx.growth = findStrPos(cnap.reacID,'BIOMASS_Ec_iJO1366_core_53p95M');
idx.atpm   = findStrPos(cnap.reacID,'ATPM');
idx.prod   = findStrPos(cnap.reacID,'EX_etoh_e');
cnap = CNAsetGenericReactionData(cnap,idx.prod,'grRules','a1 and a2 or a3'); % add genes to make reaction KI-able
idx.acEx   = findStrPos(cnap.reacID,'EX_ac_e');
idx.acUp   = findStrPos(cnap.reacID,'EX_ac_e_2');
idx.glcEx  = findStrPos(cnap.reacID,'EX_glc__D_e');
idx.glcUp  = findStrPos(cnap.reacID,'EX_glc__D_e_2');
idx.glycEx = findStrPos(cnap.reacID,'EX_glyc_e');
idx.glycUp = findStrPos(cnap.reacID,'EX_glyc_e_2');
idx.succEx = findStrPos(cnap.reacID,'EX_succ_e');
idx.succUp = findStrPos(cnap.reacID,'EX_succ_e_2');
idx.o2     = findStrPos(cnap.reacID,'EX_o2_e');

%% Set in- and out fluxes
cnap.reacMin(idx.acUp) = -30;
cnap.reacMax(idx.acUp) = 0;
cnap.reacMin(idx.acEx) = 0;
cnap.reacMax(idx.acEx) = 1000;

cnap.reacMin(idx.succUp) = 0; % Succ uptake off
cnap.reacMax(idx.succUp) = 0;
cnap.reacMin(idx.succEx) = 0;
cnap.reacMax(idx.succEx) = 1000;

cnap.reacMin(idx.glycUp) = -20;
cnap.reacMax(idx.glycUp) = 0;
cnap.reacMin(idx.glycEx) = 0;
cnap.reacMax(idx.glycEx) = 0;

cnap.reacMin(idx.glcUp)  = -10;
cnap.reacMax(idx.glcUp)  = 0;
cnap.reacMin(idx.glcEx)  = 0;
cnap.reacMax(idx.glcEx)  = 0;
cnapBU = cnap;

%% Target and desired regions
minYield = 0.3; % change to 0.4 or 0.3
Cp = 2;
Cs = 6;
Ca = 2;
Cg = 3;
T0 = full(sparse( [1          1               1              1            2          ], ...
                  [idx.prod   idx.glcUp       idx.acUp       idx.glycUp   idx.atpm   ], ...
                  [Cp         minYield*Cs     minYield*Ca    minYield*Cg  -1         ],2,cnap.numr));
t0 =  [  0  ; -3.15];
% for scenario with AcEx as notknockable
T1 = full(sparse( [1          1              1              1            2           3          ], ...
                  [idx.prod   idx.glcUp      idx.acUp       idx.glycUp   idx.acEx    idx.atpm   ], ...
                  [Cp         minYield*Cs    minYield*Ca    minYield*Cg  1           -1         ],3,cnap.numr));
t1 =  [  0  ; 0 ; -3.15];
T2 = full(sparse( [1          1                             1            2           ], ...
                  [idx.prod   idx.glcUp                     idx.glycUp   idx.atpm    ], ...
                  [Cp         minYield*Cs                   minYield*Cg  -1          ],2,cnap.numr));
t2 =  [  0  ; -3.15];
D1 = full(sparse( [1          2 ], ...
                 [idx.growth idx.atpm], ...
                 [-1         -1],2,cnap.numr));
d1 = [-0.05 ; -3.15];
D2 = full(sparse( 1   , ...
                 [idx.atpm], ...
                 -1   ,1,cnap.numr));
d2 = -8;
% yet unused


%% Knockable/Notknockable definiton
% reaction-based
[cnap,enzymes,genes,gr_rules] = CNAgenerateGERassociation(cnap);

notknockable = findStrPos(cnap.reacNotes,'.*(port|diffusion|permease|efflux|ABC).*','regex');
noGene = setdiff(1:cnap.numr,[gr_rules(:).reaction]);
spont = [gr_rules(cellfun(@(x) ismember(find(strcmp('s0001',genes)),x) , {gr_rules.genes})).reaction];
% Axel
notknockable = union(notknockable,noGene);
% new
koCost                                  = ones(cnap.numr,1);
koCost(notknockable)                    = nan;
koCost(idx.o2)                          = 1; % keep oxygen knockable
kiCost                                  = nan(cnap.numr,1);
kiCost([idx.acUp idx.glcUp idx.glycUp]) = 0.1;

% gene-based
% Axel: notKnock
% new
rkoCost                         = nan(cnap.numr,1);
rkoCost(idx.o2)                 = 1;
rkiCost                         = nan(cnap.numr,1);
rkiCost([idx.acUp idx.glcUp idx.glycUp]) = 0.1;
gkoCost                         = ones(length(genes),1);
gkoCost(strcmp(genes,'s0001'))  = nan;
gkiCost                         = nan(length(genes),1);
gkiCost(findStrPos(genes,{'a1' 'a2' 'a3'})) = 0.15; % KI of EX etoh

% computation parameters
max_num_interv = 5;
time_limit = inf;
max_solutions = inf;

%% MCS computation - reaction based
% classical
cnap.reacMin(idx.acUp)   = -30;
cnap.reacMin(idx.glycUp) =   0;
cnap.reacMin(idx.glcUp)  =   0;
mcs_ac          = CNAregMCSEnumerator(cnap,T0,t0,D1,d1,setdiff(notknockable,idx.o2),max_solutions, max_num_interv);
if ~isempty(mcs_ac), mcs_ac(:,idx.acUp) = 1; end
cnap.reacMin(idx.acUp)   =   0;
cnap.reacMin(idx.glycUp) = -15;
mcs_glyc        = CNAregMCSEnumerator(cnap,T0,t0,D1,d1,setdiff(notknockable,idx.o2),max_solutions, max_num_interv);
if ~isempty(mcs_glyc), mcs_glyc(:,idx.glycUp) = 1; end
cnap.reacMin(idx.glycUp) =   0;
cnap.reacMin(idx.glcUp)  = -10;
mcs_glc         = CNAregMCSEnumerator(cnap,T0,t0,D1,d1,setdiff(notknockable,idx.o2),max_solutions, max_num_interv);
if ~isempty(mcs_glc), mcs_glc(:,idx.glcUp) = 1; end
cnap.reacMin(idx.acUp)   = -30;
cnap.reacMin(idx.glycUp) = -15;
cnap.reacMin(idx.glcUp)  = -10;
% unite all mcs (minimal):
mcs_classic = unique([mcs_ac;mcs_glyc;mcs_glc],'rows');
[~,minimal] = mcs_isContained(mcs_classic',mcs_classic');
mcs_classic = mcs_classic(~any(tril(minimal) == 1,2),:);

%% new
[mcs_new, status, obj] = ...
   CNAmultiImposedConstrISfinder(cnap,{T1,T2},{t1,t2}, {D1},{d1},koCost,kiCost,max_num_interv+0.3, time_limit,max_solutions,1,0,1);

%% Categorize into single knock-ins or knock-in combinations
mcs_new_ac   =       find(mcs_new(idx.acUp,:)==1      & isnan(mcs_new(idx.glcUp,:))    & isnan(mcs_new(idx.glycUp,:)));
mcs_new_glc  =       find(isnan(mcs_new(idx.acUp,:))  & mcs_new(idx.glcUp,:)==1        & isnan(mcs_new(idx.glycUp,:)));
mcs_new_glyc =       find(isnan(mcs_new(idx.acUp,:))  & isnan(mcs_new(idx.glcUp,:))    & mcs_new(idx.glycUp,:)==1 );
mcs_new_glyc_ac =    find(mcs_new(idx.acUp,:)==1      & isnan(mcs_new(idx.glcUp,:))    & mcs_new(idx.glycUp,:)==1 );
mcs_new_glc_ac =     find(mcs_new(idx.acUp,:)==1      & mcs_new(idx.glcUp,:)==1        & isnan(mcs_new(idx.glycUp,:)) );
mcs_new_glyc_glc =   find(isnan(mcs_new(idx.acUp,:))  & mcs_new(idx.glcUp,:)==1        & mcs_new(idx.glycUp,:)==1 );

mcs_new_ac   = mcs_new(:,mcs_new_ac) ~= 0 & ~isnan(mcs_new(:,mcs_new_ac));
mcs_new_glc  = mcs_new(:,mcs_new_glc) ~= 0 & ~isnan(mcs_new(:,mcs_new_glc));
mcs_new_glyc = mcs_new(:,mcs_new_glyc) ~= 0 & ~isnan(mcs_new(:,mcs_new_glyc));
mcs_new_glyc_ac = mcs_new(:,mcs_new_glyc_ac) ~= 0 & ~isnan(mcs_new(:,mcs_new_glyc_ac));
mcs_new_glc_ac = mcs_new(:,mcs_new_glc_ac) ~= 0 & ~isnan(mcs_new(:,mcs_new_glc_ac));
mcs_new_glyc_glc = mcs_new(:,mcs_new_glyc_glc) ~= 0 & ~isnan(mcs_new(:,mcs_new_glyc_glc));

%% MCS computation - gene based
% Axel
% classical
cnap.reacMin(idx.acUp)   = -30;
cnap.reacMin(idx.glycUp) =   0;
cnap.reacMin(idx.glcUp)  =   0;
[gmcs_ac, gene_idx_ac]   = CNAgeneMCSEnumerator(cnap,T0, t0, D1, d1, notknockable,max_solutions, max_num_interv, [], inf, 1000, enzymes);
cnap.reacMin(idx.acUp)   =   0;
cnap.reacMin(idx.glycUp) = -15;
[gmcs_glyc, gene_idx_glyc]=CNAgeneMCSEnumerator(cnap,T0, t0, D1, d1, notknockable,max_solutions, max_num_interv, [], inf, 1000, enzymes);
cnap.reacMin(idx.glycUp) =   0;
cnap.reacMin(idx.glcUp)  = -10;
[gmcs_glc, gene_idx_glc] = CNAgeneMCSEnumerator(cnap,T0, t0, D1, d1, notknockable,max_solutions, max_num_interv, [], inf, 1000, enzymes);
cnap.reacMin(idx.acUp)   = -30;
cnap.reacMin(idx.glycUp) = -15;
cnap.reacMin(idx.glcUp)  = -10;
% unite all mcs (minimal):
gmcs_classic = unique([gmcs_ac;gmcs_glyc;gmcs_glc],'rows');
[~,minimal] = mcs_isContained(gmcs_classic',gmcs_classic');
gmcs_classic = gmcs_classic(~any(tril(minimal) == 1,2),:);

% new
[gmcs1, gcnap1, cmp_gmcs1, cmp_gcnap1, gmcs_idx1] = CNAgeneMultiImposedMCSFinder(cnap, {T1,T2} , {t1,t2} , {D1} , {d1} ,...
                                                rkoCost,rkiCost, ... koCost, kiCost
                                                max_num_interv+0.5,time_limit,max_solutions,...
                                                1,0, ... use_compression,use_bigM,
                                                1,gkoCost,gkiCost,[],1); % enum_method, gkoCost, gkiCost, enzymes, verbose
snd();
return;