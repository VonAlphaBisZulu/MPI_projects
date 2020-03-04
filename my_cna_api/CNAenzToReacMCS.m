function mcs = CNAenzToReacMCS( cnap, mcs, enzymes )
if nargin == 2
    enzymes = cnap.enzymes;
end
enzReac = {enzymes(:).reactions};
numGenes = max([enzymes(:).genes]);
for i = 1:size(enzReac,2)
    assoMat = [assoMat, iv(numGenes,[enzymes(i).genes])];
end
assoMat = assoMat';
mcs = expand_mcs(mcs, assoMat);
end
