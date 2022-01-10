%% Initialize and load models
if ~exist('cnan','var')
    startcna(1)
end
load(which('e_coli_core.mat'));
ecc = CNAcobra2cna(e_coli_core,0);
ecc2 = CNAsbmlModel2MFNetwork(which('ecc2.xml'));
ecc2_reacs = cellstr(ecc2.reacID(:,3:end));
ecc2_reacs = strrep(ecc2_reacs,'LPAREN_','');
ecc2_reacs = strrep(ecc2_reacs,'_RPAREN_','');
ecc2_reacs = strrep(ecc2_reacs,'DASH','');
ecc2_reacs = strrep(ecc2_reacs,'PPKr','PPK');
ecc2_reacs = strrep(ecc2_reacs,'SULRi','SURL');
ecc2_reacs = strrep(ecc2_reacs,'SUCCt3pp','SUCCt2_3pp');

load(which('iML1515.mat'));
load(which('iML1515geneNames.mat'));
cnap = CNAcobra2cna(iML1515,0);
iML_reacs = cellstr(cnap.reacID);

% indices
idx.glc = find(ismember(cellstr(cnap.reacID),'EX_glc__D_e'));
idx.glyc = find(ismember(cellstr(cnap.reacID),'EX_glyc_e'));
idx.ac = find(ismember(cellstr(cnap.reacID),'EX_ac_e'));
idx.fhl = find(ismember(cellstr(cnap.reacID),'FHL'));

%
reacs_keep = {};
reacs_remove = {};

% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
    ecc.reacNotes = strrep(ecc.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
gpr_rules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
[~,~,genes_ecc] = CNAgenerateGPRrules(ecc);

% block non-standard fermentation products
cnap = block_non_standard_products(cnap);
cnap.reacMin(idx.glc) = -10;
% beta-oxidation irreversible and POR5 blocked
cnap.reacMin(~cellfun(@isempty,regexp(gpr_rules,'fadE','match'))) = 0;
cnap.reacMax(ismember(cnap.reacID,{'POR5'})) = 0;
cnap.reacMin(ismember(cnap.reacID,{'POR5'})) = 0;

%% remove blocked reactions (but keep in mind that at some point glyc or ac uptake may be used)
cnap.reacMin(idx.glyc) = -10;
cnap.reacMin(idx.ac) = -10;
cnap.reacMax(idx.fhl) = 1000;
[lb, ub] = CNAfluxVariability(cnap,[],[],-1,1:cnap.numr,[],[],0);
cnap.reacMin(idx.glyc) = 0;
cnap.reacMin(idx.ac) = 0;
cnap.reacMax(idx.fhl) = 0;
r_remove = abs(lb)<1e-9 & abs(ub)<1e-9 ;
reacs_remove = unique([reacs_remove; cellstr(cnap.reacID(r_remove,:))]);

% reduce model for the first time, get new indices
cnap = CNAdeleteReaction(cnap,find(ismember(cellstr(cnap.reacID),reacs_remove)));
cnap = CNAdeleteSpecies(cnap,find(~any(cnap.stoichMat,2)),0);
gpr_rules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
idx.glc = find(ismember(cellstr(cnap.reacID),'EX_glc__D_e'));
idx.glyc = find(ismember(cellstr(cnap.reacID),'EX_glyc_e'));
idx.ac = find(ismember(cellstr(cnap.reacID),'EX_ac_e'));
idx.bm = find(~cellfun(@isempty,regexp(cellstr(cnap.reacID),'BIOMASS_Ec.*core.*','match')));
idx.etoh = find(ismember(cellstr(cnap.reacID),'EX_etoh_e'));
idx.lac = find(ismember(cellstr(cnap.reacID),'EX_lac__D_e'));
idx.succ = find(ismember(cellstr(cnap.reacID),'EX_succ_e'));
idx.o2 = find(ismember(cellstr(cnap.reacID),'EX_o2_e'));
idx.fhl = find(ismember(cellstr(cnap.reacID),'FHL'));

%% keep all reactions that are catalyzed by genes from e_coli_core
% r_keep = false(cnap.numr,1);
% for i = 1:numel(genes_ecc)
%     r_keep = r_keep | ~cellfun(@isempty,regexp(gpr_rules,genes_ecc(i)));
% end
% reacs_keep = [reacs_keep; cellstr(cnap.reacID(r_keep,:))];

%% keep all reactions that occur in ECC2 and iML1515
reacs_keep = [reacs_keep; iML_reacs(ismember(iML_reacs,ecc2_reacs))];

%% keep reactions for maximal aerobic growth
fv = CNAoptimizeFlux(cnap,[], [], -1, -1, 1); % pFBA
% V = zeros(1,cnap.numr);
% V(idx.bm) = -1;
% v = -fv(idx.bm);
% [lb, ub] = CNAfluxVariability(cnap,[],[],-1,1:cnap.numr,V,v,0);
% keep all reactions that are essential for maximal growth (aerobic)
%    % if lb and ub have same sign, reactions are essential for growth
% keep all reactions from max-growth pFBA
% r_keep = (lb.*ub)>1e-8 | abs(fv)>1e-8;
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% keep reactions for maximal anaerobic growth
cnap_anaerobic = cnap;
cnap_anaerobic.reacMin(idx.o2) = 0;
cnap_anaerobic.reacMax(idx.o2) = 0;
fv = CNAoptimizeFlux(cnap_anaerobic,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% keep reactions for maximal anaerobic ethanol production
cnap_anaerobic.objFunc(:) = 0;
cnap_anaerobic.objFunc(idx.etoh) = -1;
fv = CNAoptimizeFlux(cnap_anaerobic,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% keep reactions for maximal anaerobic lactate production
cnap_anaerobic.objFunc(:) = 0;
cnap_anaerobic.objFunc(idx.lac) = -1;
fv = CNAoptimizeFlux(cnap_anaerobic,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% keep reactions for maximal anaerobic acetate production
cnap_anaerobic.objFunc(:) = 0;
cnap_anaerobic.objFunc(idx.ac) = -1;
fv = CNAoptimizeFlux(cnap_anaerobic,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);
% keep also POX
reacs_keep = [reacs_keep; {'POX'}];

%% keep reactions for maximal anaerobic succinate production
cnap_anaerobic.objFunc(:) = 0;
cnap_anaerobic.objFunc(idx.succ) = -1;
fv = CNAoptimizeFlux(cnap_anaerobic,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);
% include different succinate transporters
reacs_keep = [reacs_keep; {'SUCCt1pp';'SUCCt2_2pp'}];

%% keep reactions for glycerate uptake
cnap_glyc = cnap;
cnap_glyc.reacMin(idx.glyc) = -10;
cnap_glyc.reacMax(idx.glyc) = -10;
cnap_glyc.reacMin(idx.glc) = 0;
cnap_glyc.reacMax(idx.glc) = 0;
% [lb, ub] = CNAfluxVariability(cnap_glyc,[],[],-1,1:cnap.numr,[],[],0);
% r_keep = (lb.*ub)>1e-8;
fv = CNAoptimizeFlux(cnap_glyc,[], [], -1, -1, 1); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);
% % keep glycerol dehydrogenase (NADH), glycerol kinase,
% % glycerol-3-phosphate dehydrogenase (ubiquinone-8, NADPH)
% reacs_keep = unique([reacs_keep; {'GLYCDx';'GLYK';'G3PD2';'G3PD5'}]);

% %% keep reactions for FHL
% cnap_fhl = cnap;
% cnap_fhl.reacMin(idx.o2) = 0;
% cnap_fhl.reacMin(idx.fhl) = 1;
% cnap_fhl.reacMax(idx.fhl) = 1;
% cnap_fhl.objFunc(:) = 0;
% % [lb, ub] = CNAfluxVariability(cnap_glyc,[],[],-1,1:cnap.numr,[],[],0);
% % r_keep = (lb.*ub)>1e-8;
% fv = CNAoptimizeFlux(cnap_fhl,[], [], -1, -1, 1); % pFBA
% r_keep = abs(fv)>1e-8;
% reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

% %% keep ASPT pathway
% idx.aspt = find(ismember(cellstr(cnap.reacID),'ASPT'));
% cnap.objFunc(:) = 0;
% cnap.objFunc(idx.aspt) = -1;
% fv = CNAoptimizeFlux(cnap,[], [], -1, -1, 1); % pFBA
% r_keep = abs(fv)>1e-8;
% reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

% %% keep mththf pathway
% idx.aspt = find(ismember(cellstr(cnap.reacID),'ASPT'));
% cnap.objFunc(:) = 0;
% cnap.objFunc(idx.succ) = -1;
% fv = CNAoptimizeFlux(cnap,[], [], -1, -1, 1); % pFBA
% r_keep = abs(fv)>1e-8;
% reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% fab-gene-Reactions
% Simplification: Find reactions that are strictly positive and attain high flux values. 
% Demand that they run in forward direction
fab_reacs = find(~cellfun(@isempty,regexp(gpr_rules,'fab','match')));
[lb, ub] = CNAfluxVariability(cnap,[],[],-1,fab_reacs,[],[],0);
pos_fab_reacs = fab_reacs(abs(lb)<=1e-9 & ub>=0.1);
V = sparse(1:numel(pos_fab_reacs),pos_fab_reacs,-1,numel(pos_fab_reacs),cnap.numr);
v = -0.01*ones(numel(pos_fab_reacs),1);
fv = CNAoptimizeFlux(cnap,[], [], -1, -1, 1,V,v); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);
% [lb, ub] = CNAfluxVariability(cnap,[],[],-1,1:cnap.numr,V,v,0);
% r_keep = (lb.*ub)>1e-8;
% reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% fad-gene-Reactions
% Simplification: Use pFBA and enforce positive sum of fluxes to make sure one degradation pathway is feasible
fad_reacs = find(~cellfun(@isempty,regexp(gpr_rules,'fad','match')));
[lb, ub] = CNAfluxVariability(cnap,[],[],-1,fad_reacs,[],[],0);
pos_fad_reacs = fad_reacs(abs(lb)<=1e-9 & ub>=0.1);
V = sparse(1:numel(pos_fad_reacs),pos_fad_reacs,-1,numel(pos_fad_reacs),cnap.numr);
v = -0.01*ones(numel(pos_fad_reacs),1);
fv = CNAoptimizeFlux(cnap,[], [], -1, -1, 1,V,v); % pFBA
r_keep = abs(fv)>1e-8;
reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);
% [lb, ub] = CNAfluxVariability(cnap,[],[],-1,1:cnap.numr,V,v,0);
% r_keep = (lb.*ub)>1e-8;
% reacs_keep = unique([reacs_keep; cellstr(cnap.reacID(r_keep,:))]);

%% Make preliminary core get gene-list and make sure all reactions with the same genes are included in the final core model
cnap_core_prelim = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),reacs_keep)));
[~,~,genes_prelim] = CNAgenerateGPRrules(cnap_core_prelim);
r_keep = false(cnap.numr,1);
for i = 1:numel(genes_prelim)
%     disp(sum(~cellfun(@isempty,regexp(gpr_rules,genes_prelim(i)))));
%     disp(genes_prelim(i))
    r_keep = r_keep | ~cellfun(@isempty,regexp(gpr_rules,genes_prelim(i)));
end
return
reacs_keep = [reacs_keep; cellstr(cnap.reacID(r_keep,:))];

[lb, ub] = CNAfluxVariability(cnap_core_prelim,[],[],-1,1:cnap_core_prelim.numr,[],[],0);

return;
%% Make core model and check feasibility of all reactions with FVA
cnap_core = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),reacs_keep)));
[lb, ub] = CNAfluxVariability(cnap_core,[],[],-1,1:cnap_core.numr,[],[],0);

% cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
% cnap = CNAdeleteSpecies(cnap,find(~any(cnap.stoichMat,2)),0);