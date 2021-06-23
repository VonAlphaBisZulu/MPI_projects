if ~exist('cnan','var')
    startcna(1)
end
% load CNA network
CNAloadNetwork({'network_dirs/iJOcore';1});
% load product pathway
[product_rID,species,reactions] = load_pathway(12); % 12: bisabolene
for spec = species
    iJOcore = CNAaddSpeciesMFN(iJOcore,spec.spec_id,0,spec.spec_name);
    iJOcore = CNAsetGenericSpeciesData(iJOcore,iJOcore.nums,'fbc_chemicalFormula',char(spec.fbc_chemicalFormula),'fbc_charge',double(spec.fbc_charge));
end
for reac = reactions
    iJOcore = CNAaddReactionMFN(iJOcore, reac.reac_id, reac.equation, reac.lb, reac.ub,0,nan,nan,'',0,0,0,0);
    iJOcore = CNAsetGenericReactionData(iJOcore,iJOcore.numr,'geneProductAssociation',char(reac.fbc_geneProductAssociation));
end
iJOcore.reacBoxes(iJOcore.reacBoxes(:,5) == 0,5) = -1;
iJOcore = CNAgenerateMap(iJOcore,1,0);
return
load(which('EX_bsb_e_iJO.mat'));
% switch rows for correct network mapping
netmap = findStrPos(cellstr(iJOcore.reacID),cellstr(cnap.reacID));
rmcs(netmap,:) = rmcs(1:numel(netmap),:);
CNAloadMCSinGUI(iJOcore,rmcs);
