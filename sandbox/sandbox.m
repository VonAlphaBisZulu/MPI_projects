cnap = CNAloadNetwork({'../network_dirs/yield_net';1},1,1,0);
r_P = 5;
r_S = 6;
r_BM = 7;
modules = {struct};
modules{1}.type = 'yield_w_constr';
modules{1}.sense = 'target';
% r*e/r*f <= y
modules{1}.e = zeros(1,cnap.numr);
modules{1}.e(r_P) = 1; % e numerator
modules{1}.f = zeros(1,cnap.numr);
modules{1}.f(r_BM) = 1; % f denominator
modules{1}.y = 0;
modules{1}.V = [];
modules{1}.v = [];

modules{2}.type = 'lin_constraints';
modules{2}.sense = 'desired';
modules{2}.V = zeros(1,cnap.numr);
modules{2}.V(r_BM) = -1; % f denominator
modules{2}.v = -0.5;

obj = MCS_cplex(cnap,modules);
[mcs,status] = obj.findSmallestMCS(maxSolutions,inf);
CNAloadMCSinGUI(yield_net,-mcs);