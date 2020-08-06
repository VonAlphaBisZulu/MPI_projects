function mcs = CNAgeneToReacMCS( cnap, enzymes, gmcs, gidx)
if isempty(enzymes)
    enzymes = cnap.enzymes;
end

% expand gMCS
numGenes = max([enzymes(:).genes]); % get total number of Genes in enzyme structure
ggMCSMap = sparse(1:length(gidx),gidx,1,length(gidx),numGenes);
gmcs = gmcs*ggMCSMap;

% mapping to eMCS
geneEnz = {enzymes(:).genes};
geMap = [];
for i = 1:length(geneEnz)
    geMap(geneEnz{i},i) = 1;
end
emcs = gmcs*geMap;

% mapping to reaction MCS
enzReac = {enzymes(:).reactions};
erMap = [];
for i = 1:length(enzReac)
    erMap(i,enzReac{i}) = 1;
end

% translate to reaction cuts and only count those where every possible enzyme
% is cut
numEnzPerReacMat = repmat(sum(erMap),size(emcs,1),1);
numNonZeroEnzPerReacMat = repmat(sum(erMap)>0,size(emcs,1),1);

mcs = ((emcs*erMap) == numEnzPerReacMat) & numNonZeroEnzPerReacMat;
end
