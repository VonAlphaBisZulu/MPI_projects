if ~exist('cnan','var')
    startcna(1)
end
% load Bigg network
load(which('iML1515.mat'));
cnap = CNAcobra2cna(iML1515,0);
load(which('core.mat'));
cnap = CNAdeleteReaction(cnap,find(~ismember(cellstr(cnap.reacID),core_reacs)));
cnap = CNAdeleteSpecies(cnap,find(~ismember(cellstr(cnap.specID),core_specs)),0);
clear iML1515;
% load CNA network
CNAloadNetwork({'network_dirs/iMLcore';1});
load(which('iML1515geneNames.mat'));
iMLcore.reacMin(ismember(iMLcore.reacID,{'EX_glc__D_e'})) = -10;
% remove non-core reactions
iMLcore = CNAdeleteReaction(iMLcore,find(~ismember(cellstr(iMLcore.reacID),core_reacs)));
iMLcore = CNAdeleteSpecies(iMLcore,find(~ismember(cellstr(iMLcore.specID),core_specs)),0);
% replace gene numbers with names
for i = 1:length(ecoliGeneNames)
    iMLcore.reacNotes = strrep(iMLcore.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
end
% load product pathway
[product_rID,species,reactions] = load_pathway(12); % 12: bisabolene
for spec = species
    iMLcore = CNAaddSpeciesMFN(iMLcore,spec.spec_id,0,spec.spec_name);
    iMLcore = CNAsetGenericSpeciesData(iMLcore,iMLcore.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    iMLcore = CNAaddReactionMFN(iMLcore, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    iMLcore = CNAsetGenericReactionData(iMLcore,iMLcore.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end
iMLcore = CNAsetGenericReactionData_with_array(iMLcore,'on_generated_map',num2cell(zeros(iMLcore.numr,1)));
iMLcore.reacBoxes(iMLcore.reacBoxes(:,5) == 0,5) = -1;
iMLcore = CNAgenerateMap(iMLcore,1,0);
load(which('EX_bsb_e_core.mat'));
% switch rows for correct network mapping
netmap = findStrPos(cellstr(iMLcore.reacID),cellstr(cnap.reacID));% map networks
rmcs_tot(netmap,:) = rmcs_tot(1:numel(netmap),:);
CNAloadMCSinGUI(iMLcore,rmcs_tot);
