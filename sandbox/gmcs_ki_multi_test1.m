% 02.07.2019
% Philipp
% Testing the use of multiple substrates and multiple target and desired regions in iJO1366 with gene MCS:
%
%% Preparation
cnap = CNAloadNetwork({'iJO1366';1},1,1);

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

% Set upper and lower bounds for sources and sinks
% find exchange reactions (reactions that have only one enry, and that entry is -1)
ex_reacs = find(sum((sum(abs(cnap.stoichMat))==1).*(cnap.stoichMat==-1)));

specsWithCarbon = regexp(cellstr(cnap.specNotes), '\[.*C([A-Z]|\d).*]', 'match');
specsWithCarbon = find(~cellfun(@isempty,specsWithCarbon));
reacsWCarbon  = cnap.reacID(ex_reacs( ismember(ex_reacs,find(sum(cnap.stoichMat(specsWithCarbon,:))))),:);

% RMin

cnap.reacMin(findStrPos(cnap.reacID,reacsWCarbon)) = 0;
cnap.reacMin(findStrPos(cnap.reacID,'EX_glc__D_e')) = -10;
cnap.reacMin(findStrPos(cnap.reacID,'EX_co2_e')) = -1000;

% RMax
cnap.reacMax(ex_reacs) = 0;
cnap.reacMax(findStrPos(cnap.reacID,'BIOMASS.*core','regex')) = 1000;
cnap.reacMax(findStrPos(cnap.reacID,'BIOMASS.*WT','regex')) = 0;
cnap.reacMax(findStrPos(cnap.reacID,{   'DM_4crsol_c'...
                                        'DM_5drib_c'...
                                        'DM_aacald_c'...
                                        'DM_amob_c'...
                                        'DM_mththf_c'...
                                        'DM_oxam_c'...
                                        'EX_ac_e'...
                                        'EX_co2_e'...
                                        'EX_etoh_e'...
                                        'EX_for_e'...
                                        'EX_glyc_e'...
                                        'EX_glyc__R_e'...
                                        'EX_h2_e'...
                                        'EX_h2o2_e'...
                                        'EX_h2o_e'...
                                        'EX_h_e'...
                                        'EX_lac__D_e'...
                                        'EX_meoh_e'...
                                        'EX_o2_e'...
                                        'EX_succ_e'...
                                        'EX_tungs_e'})) = 1000;
                                    
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
%% T0, T1 and T2
minYield = 0.3; % change to 0.4 or 0.3
Cp = 2;
Cs = 6;
Ca = 2;
Cg = 3;
T0 = full(sparse( [1          1                             1            2          ], ...
                  [idx.prod   idx.glcUp                     idx.glycUp   idx.atpm   ], ...
                  [Cp         minYield*Cs                   minYield*Cg  -1         ],2,cnap.numr));
t0 =  [  0  ; -3.15];
% for scenario with AcEx as notknockable
T1 = full(sparse( [1          1              1              1            2           3          ], ...
                  [idx.prod   idx.glcUp      idx.acUp       idx.glycUp   idx.acEx    idx.atpm   ], ...
                  [Cp         minYield*Cs    minYield*Ca    minYield*Cg  1           -1         ],3,cnap.numr));
t1 =  [  0  ; 0 ; -3.15];
T2 = full(sparse( [1          1                             1            2          ], ...
                  [idx.prod   idx.glcUp                     idx.glycUp   idx.atpm   ], ...
                  [Cp         minYield*Cs                   minYield*Cg  -1         ],2,cnap.numr));
t2 =  [  0  ; -3.15];

%% D
% Desired: Growth 0.1 at ATPM = 1.5
D1 = full(sparse( [1          2 ], ...
                 [idx.growth idx.atpm], ...
                 [-1         -1],2,cnap.numr));
d1 = [-0.05 ; -1.5];

D2 = full(sparse( 1   , ...
                 [idx.atpm], ...
                 -1   ,1,cnap.numr));
d2 = -8;

[cnap,enzymes,genes,gr_rules] = CNAgenerateGERassociation(cnap);
notknockable = full(sparse(findStrPos(cnap.reacNotes,'.*(port|diffusion|permease|efflux|ABC).*','regex'),1,1,cnap.numr,1));
noGene = cellfun(@isempty,CNAgetGenericReactionData_as_array(cnap,'grRules'));
spont = cellfun(@(x) strcmp('s0001',x),CNAgetGenericReactionData_as_array(cnap,'grRules'));
notknockable = find(spont | noGene | notknockable);

% Axel
notknockable = union(notknockable,find(noGene));
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

max_num_interv = 6;
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
                                                1,gkoCost,gkiCost,[],0); % enum_method, gkoCost, gkiCost, enzymes, verbose

return
koCost(idx.acEx) = nan; % make Acetate export notknockable again
is=is1;
% % AcEx is notknockable, but there are two target regions
% [is, status, obj] = ...
%    CNAmultiImposedConstrISfinder(cnap,{T1,T2},{t1,t2},  {D1,D2},{d1,d2},koCost,kiCost,max_num_interv,   time_limit,max_solutions,1,0,1);
% 
% is1_bu = is1;
is_bu = is;

%% compare the IS with two target regions and one (with the option of knocking-out AcEx)
is1 = is1(:,~(isnan(is1(7,:)) & (is1(3,:)==-1))); % remove all IS with Ac-KO but not Ac-KI
is1 = is1(:,sum(~isnan(is1) & (is1 ~= 0),1)<=max_num_interv); % remove all IS with too many interventions
is1(idx.acEx,:) = 0;
num_total_is = size(unique([~isnan(is1) & is1 ~= 0 , ~isnan(is) & is ~= 0]','rows')',2);
if num_total_is == size(is,2), disp('is-sets are identical'); else, disp('is-sets are not identical'); end
is1 = is1_bu;

%% Categorize into single knock-ins or knock-in combinations
is_ac   = find(is(7,:)==1      & isnan(is(34,:))  & isnan(is(36,:))) ;
is_ac   = is(:,is_ac) == -1;
is_glc  = find(isnan(is(7,:))  & is(34,:)==1      & isnan(is(36,:))) ;
is_glc  = is(:,is_glc) == -1;
is_glyc = find(isnan(is(7,:))  & isnan(is(34,:))  & is(36,:)==1    );
is_glyc = is(:,is_glyc) == -1;
is_glyc_ac = find(is(7,:)==1      & isnan(is(34,:))  & is(36,:)==1    );
is_glyc_ac = is(:,is_glyc_ac) == -1;
is_glc_ac = find(is(7,:)==1       & is(34,:)==1   & isnan(is(36,:))   );
is_glc_ac = is(:,is_glc_ac) == -1;
is_glyc_glc = find(isnan(is(7,:))      & is(34,:)==1  & is(36,:)==1    );
is_glyc_glc = is(:,is_glyc_glc) == -1;
is_glyc_glc_ac = find(is(7,:)==1      & is(34,:)==1   & is(36,:)==1    );
is_glyc_glc_ac = is(:,is_glyc_glc_ac) == -1;

num_intv = sum((is==-1 | is==1),1);
%% Traditional MCS computation
T3 = full(sparse( [1          1                                          2           3           4            5         ], ...
                  [idx.prod   idx.glcUp                                  idx.glcEx   idx.acUp    idx.glycUp   idx.atpm  ], ...
                  [Cp         minYield*Cs                                1           1           1            -1        ],5,cnap.numr));
t3 =  [  0  ; 0 ; 0 ; 0 ; -3.15];
T4 = full(sparse( [1                         1                           2           3           4            5         ], ...
                  [idx.prod                  idx.acUp                    idx.glcUp   idx.acEx    idx.glycUp   idx.atpm  ], ...
                  [Cp                        minYield*Ca                 1           1           1            -1        ],5,cnap.numr));
t4 =  [  0  ; 0 ; 0 ; 0 ; -3.15];
T5 = full(sparse( [1                                        1            2           3           4            5         ], ...
                  [idx.prod                                 idx.glycUp   idx.glcUp   idx.acUp    idx.glycEx   idx.atpm  ], ...
                  [Cp                                       minYield*Cg  1           1           1            -1        ],5,cnap.numr));
t5 =  [  0  ; 0 ; 0 ; 0 ; -3.15];

% make sure that uptake-reactions cannot be knocked in / out
notknockable = isnan(koCost); 

%% Test KO Strategies for Glc only
cnap = cnapBU;
cnap.reacMin(idx.acUp)   = 0;
cnap.reacMax(idx.acUp)   = 0;
cnap.reacMin(idx.glycUp) = 0;
cnap.reacMax(idx.glycUp) = 0;
cnap.reacMin(idx.glcUp)  = -10;
cnap.reacMax(idx.glcUp)  = 0;

mcs_glc = CNAregMCSEnumerator(cnap,T3,t3,D1,d1,find(notknockable),max_solutions,max_num_interv-1,[],1);
mcs_glc = logical(mcs_glc');
a = 0;
for i = 1:size(mcs_glc,2)
    cnap = cnapBU;
    cnap.reacMin(mcs_glc(:,i)) = 0;
    cnap.reacMax(mcs_glc(:,i)) = 0;
    a(i) = CNAcheckRegionFeasibility(cnap,D2,d2);
end
mcs_glc = mcs_glc(:,logical(a));
%% Test KO Strategies for Ac only
cnap = cnapBU;
cnap.reacMin(idx.acUp)   = -30;
cnap.reacMax(idx.acUp)   = 0;
cnap.reacMin(idx.glycUp) = 0;
cnap.reacMax(idx.glycUp) = 0;
cnap.reacMin(idx.glcUp)  = 0;
cnap.reacMax(idx.glcUp)  = 0;

mcs_ac = CNAregMCSEnumerator(cnap,T4,t4,D1,d1,find(notknockable),max_solutions,max_num_interv-1,[],1);
mcs_ac = logical(mcs_ac');
a = 0;
for i = 1:size(mcs_ac,2)
    cnap = cnapBU;
    cnap.reacMin(mcs_ac(:,i)) = 0;
    cnap.reacMax(mcs_ac(:,i)) = 0;
    a(i) = CNAcheckRegionFeasibility(cnap,D2,d2);
end
mcs_ac = mcs_ac(:,logical(a));
%% Test KO Strategies for Glyc only
cnap = cnapBU;
cnap.reacMin(idx.acUp)   = 0;
cnap.reacMax(idx.acUp)   = 0;
cnap.reacMin(idx.glycUp) = -20;
cnap.reacMax(idx.glycUp) = 0;
cnap.reacMin(idx.glcUp)  = 0;
cnap.reacMax(idx.glcUp)  = 0;

mcs_glyc = CNAregMCSEnumerator(cnap,T5,t5,D1,d1,find(notknockable),max_solutions,max_num_interv-1,[],0);
mcs_glyc = logical(mcs_glyc');
a = 0;
for i = 1:size(mcs_glyc,2)
    cnap = cnapBU;
    cnap.reacMin(mcs_glyc(:,i)) = 0;
    cnap.reacMax(mcs_glyc(:,i)) = 0;
    a(i) = CNAcheckRegionFeasibility(cnap,D2,d2);
end
mcs_glyc = mcs_glyc(:,logical(a));

%% Are all single knock-in IS identical to the 'traditional' MCS?
if size(mcs_ac,2) > 0 && size(is_ac,2) > 0
    if sum(ismember(is'==-1,mcs_ac','rows')')   == size(mcs_ac,2), disp('Ac alle enthalten'); else, disp('Ac nicht alle enthalten'); end
end
if size(mcs_glyc,2) > 0 && size(mcs_glyc,2) > 0
    if sum(ismember(is'==-1,mcs_glyc','rows')') == size(mcs_glyc,2), disp('Glyc alle enthalten'); else, disp('Glyc nicht alle enthalten'); end
end
if size(mcs_glc,2) > 0 && size(is_glc,2) > 0
    if sum(ismember(is'==-1,mcs_glc','rows')')  == size(mcs_glc,2), disp('Glc alle enthalten'); else, disp('Glc nicht alle enthalten'); end
end
    
% save('./_StrainBooster/sandbox.mat','is','mcs_ac','mcs_glyc','mcs_glc');
%% validation of IS
for i = 1:size(is,2)
    cnap = cnapBU;
    cnap.reacMin(is(:,i)==-1 | isnan(is(:,i))) = 0;
    cnap.reacMax(is(:,i)==-1 | isnan(is(:,i))) = 0;

    if CNAcheckRegionFeasibility(cnap,T1,t1)
    disp(['Target feasible in mcs ' num2str(i)]); else
    disp(['Target infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckRegionFeasibility(cnap,T2,t2)
    disp(['Target feasible in mcs ' num2str(i)]); else
    disp(['Target infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckRegionFeasibility(cnap,D1,d1)
    disp(['Desired feasible in mcs ' num2str(i)]); else
    disp(['Desired infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckRegionFeasibility(cnap,D2,d2)
    disp(['Desired feasible in mcs ' num2str(i)]); else
    disp(['Desired infeasible in mcs ' num2str(i)]); 
    end

    if ~isnan(is(idx.acUp,i))
        cnap.reacMin(idx.acEx) = 0;
        cnap.reacMax(idx.acEx) = 0;
    end
    
    [Ymin,~,~,~] = CNAoptimizeYield(cnap,-iv(cnap.numr,idx.prod)',-iv(cnap.numr,idx.glcUp)');
    disp(['Ymin (glc only M/M) = ' num2str(Ymin)]); 
    [Ymin,~,~,~] = CNAoptimizeYield(cnap,-Cp*iv(cnap.numr,idx.prod)',-Cs*iv(cnap.numr,idx.glcUp)'-Ca*iv(cnap.numr,idx.acUp)');
    disp(['Ymin (total C balance) = ' num2str(Ymin) 10 '______']); 
end