% This script computes yields for octyl acetate production in the 
% iJO1366 core model with blocked reverse beta-oxidation
% Philipp - 21.10.2021

%% 1. build model
if ~exist('cnan','var')
    startcna(1)
end
% load model
load(which('iJO1366.mat'));
load(which('iJO1366geneNames.mat'));
cnap = CNAcobra2cna(iJO1366,0);
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
gpr_gs = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
% eliminate reactions with blocked flux
cnap = block_non_standard_products(cnap);
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
cnap.reacMin(~cellfun(@isempty,regexp(gpr_gs,'fadE','match'))) = 0; % block reverse beta-oxidation
[minFlux,maxFlux] = CNAfluxVariability(cnap,[],[],-1,1:cnap.numr,[],[],0);
cnap = CNAdeleteReaction(cnap,find(minFlux==0 & maxFlux==0));
cnap = CNAdeleteSpecies(cnap,find(~any(cnap.stoichMat,2)),0);
% load core reaction names
load(which('iJO_core.mat'));
gpr = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
[product_rID,species,reactions] = load_pathway(14);
core_reacs = [core_reacs;{reactions.reac_id}'];
core_reacs = [core_reacs;{'NADH18pp'; 'SULR'}];
% conserve also all fad-reactions
core_reacs = [core_reacs;cellstr(cnap.reacID(find(~cellfun(@isempty,regexp(gpr,'fad','match'))),:))];
core_specs = [core_specs;{species.spec_id}'  ];
% Protect all reactions that provide octa_c
cnap1 = cnap;
cnap1.reacMin(cnap1.reacMin>=0) = 0; % deactivate ATPM
octa_reacs = find(cnap1.stoichMat(strcmp(cellstr(cnap1.specID),'octa_c'),:));
for i = octa_reacs
    cnap1.objFunc(:) = 0;
    cnap1.objFunc(i) = -1;
    fv = CNAoptimizeFlux(cnap1,[],[],2,0);
    if fv(i)>=0.1
        A_ieq = sparse(1,i,-1,1,cnap1.numr);
        b_ieq = -fv(i);
        [minFlux,maxFlux] = CNAfluxVariability(cnap1,[],[],-1,1:cnap1.numr,A_ieq,b_ieq,0);
        octa_c_essential = find((minFlux.*maxFlux)>0);
        core_reacs = unique([core_reacs;cellstr(cnap1.reacID(octa_c_essential,:))]);
    end
end
% Add heterologous pathway
for spec = species
    cnap = CNAaddSpeciesMFN(cnap,spec.spec_id,0,spec.spec_name);
    cnap = CNAsetGenericSpeciesData(cnap,cnap.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    cnap = CNAaddReactionMFN(cnap, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    cnap = CNAsetGenericReactionData(cnap,cnap.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end
% Reduce model to core network
cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
cnap = CNAdeleteSpecies(cnap,find(~any(cnap.stoichMat,2)),0);
% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% Checking mass and and charge balances after pathway additions
check_mass_balance(cnap);
% replace bounds with inf
cnap.reacMin(cnap.reacMin==-1000) = -inf;
cnap.reacMax(cnap.reacMax== 1000) =  inf;
% lump gene names
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

%% Compute yields
substrate_rID = 'EX_glc__D_e';
biomass_rID   = regexp(cellstr(cnap.reacID),'BIOMASS_Ec.*core.*','match');
biomass_rID   = char([biomass_rID{:}]);

idx.subs = find(ismember(cellstr(cnap.reacID),substrate_rID));
idx.bm = find(ismember(cellstr(cnap.reacID),biomass_rID));
idx.prod = find(ismember(cellstr(cnap.reacID),product_rID));

gpr = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
% to block ethanol production, actually adhE, adhP, frmA, mhpF
reac_adhE = find(~cellfun(@isempty,regexp(gpr,'mhpF or adhE','match'))); % simplified 
% to block ethanol production, actually ldhA, dld
reac_ldhA = find(~cellfun(@isempty,regexp(gpr,'dld or ldhA','match')));
% oxygen exchange
reac_o2   = find(strcmp(cellstr(cnap.reacID),'EX_o2_e'));

% knockout-scenarios
% 1: no ko, aerobic
sc{1} = [];
% 2: no ko, anaerobic
sc{2} = reac_o2;
% 3: adhE, aerobic
sc{3} = reac_adhE;
% 4: adhE, anaerobic
sc{4} = [reac_o2, reac_adhE];
% 5: adhE, anaerobic
sc{5} = [reac_adhE, reac_ldhA];
% 6: adhE, anaerobic
sc{6} = [reac_o2, reac_adhE, reac_ldhA];

c = zeros(1,cnap.numr);
c(idx.prod) = 1;
d = zeros(1,cnap.numr);
d(idx.subs) = -1;
c2 = zeros(1,cnap.numr);
c2(idx.bm) = 1;

for i = 1:6
    disp('==============');
    txt = ['Scenario ' num2str(i) ', KOs: ' strjoin(cellstr(cnap.reacID(sc{i},:)),', ')];
    disp(txt);
    cnap_i = cnap;
    cnap_i.reacMin(sc{i}) = 0;
    cnap_i.reacMax(sc{i}) = 0;
    cnap_i.objFunc(idx.bm) = -1;
    fv = CNAoptimizeFlux(cnap_i);
    disp(['mue max: ' num2str(fv(idx.bm))]);
    maxyield = CNAoptimizeYield(cnap_i,c,d,[],[],2,-1);
    disp(['maxyield: ' num2str(maxyield)]);
    minyield = CNAoptimizeYield(cnap_i,-c,d,[],[],2,-1);
    disp(['minyield: ' num2str(-minyield)]);
    if i == 6
        disp('plotting..');
        [yieldspace,success,fig] = CNAplotYieldSpace(cnap_i,c2,d,c,d,20,[],[],2,0);
        title(strrep(txt,'_','\_'));
    end
end
