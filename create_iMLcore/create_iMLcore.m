load(which('iML1515.mat'));
load(which('iML1515geneNames.mat'));
load(which('e_coli_core.mat'));
cnap = CNAcobra2cna(iML1515,0);
cnap = block_non_standard_products(cnap);
cnap.reacMin(ismember(cnap.reacID,{'EX_glc__D_e'})) = -10;
cnap.path = 'iMLcore';
load(which('core.mat'));
% Reduce model to core network
cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
cnap = CNAdeleteSpecies(cnap,find(~ismember(cellstr(cnap.specID),core_specs)),0);
% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end