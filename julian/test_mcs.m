if ~exist('cnan','var')
    startcna(1)
end
% load Bigg network
load(which('iML1515.mat'));
cnap = CNAcobra2cna(iML1515,0);
clear iML1515;
% load CNA network
CNAloadNetwork({'network_dirs/iML1515';1});
load(which('iML1515geneNames.mat'));
iML1515 = block_non_standard_products(iML1515);
iML1515.reacMin(ismember(iML1515.reacID,{'EX_glc__D_e'})) = -10;
% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    iML1515.reacNotes = strrep(iML1515.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% load product pathway
[product_rID,species,reactions] = load_pathway(12); % 12: bisabolene
for spec = species
    iML1515 = CNAaddSpeciesMFN(iML1515,spec.spec_id,0,spec.spec_name);
    iML1515 = CNAsetGenericSpeciesData(iML1515,iML1515.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    iML1515 = CNAaddReactionMFN(iML1515, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    iML1515 = CNAsetGenericReactionData(iML1515,iML1515.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end
iML1515.reacBoxes(iML1515.reacBoxes(:,5) == 0,5) = -1;
iML1515 = CNAgenerateMap(iML1515,1,0);
load(which('EX_bsb_e.mat'));
CNAloadMCSinGUI(iML1515,rmcs);