load('matlab3.mat');
CNAloadNetwork({'../network_dirs/iMLcore';1});
iMLcore = CNAaddSpecsAndReacsFromFile(iMLcore,prod{:});