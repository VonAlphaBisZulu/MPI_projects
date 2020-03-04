tic
startcna(1);
model       = 'ECGS';
metabolite  = 'P18_Sorbitol_Gluconic_acid';
ECGS = CNAloadNetwork(find(strcmp(cellstr({cnan.net.path}'),model)),1,1);
ECGS = CNAaddSpecsAndReacsFromFile(ECGS,[metabolite '.xls']);
ECGS_MCS;
filename=['_StrainBooster/_My_Simulations/' model '-' metabolite '-' datestr(date,'yyyy-mm-dd') '.mat'];
save(filename,'mcs','cnap','ECGS','desiredReacs','d','D','targetReacs','t','T','genericConstrains','notknockable','maxMCSnum','maxMCSsize');
clear cnan cnap ECGS model metabolite filename mcs ECGS desiredReacs d D targetReacs t T genericConstrains notknockable maxMCSnum maxMCSsize
toc