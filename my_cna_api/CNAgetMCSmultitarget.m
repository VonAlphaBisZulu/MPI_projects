if ~exist('cnan','var')
    startcna(1)
end
cnap = CNAloadNetwork({'iJO1366';1},1,1);
filename = 'P17_Acetoin_2.xls';

cnap = CNAaddSpecsAndReacsFromFile(cnap,filename);
[T, t, D, d,notknockable,reacMin,reacMax,xor] = CNAgetMCScalcParamXls( cnap, filename);
% xor still needs to be implemented