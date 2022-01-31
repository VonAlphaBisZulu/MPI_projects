%% Initialize and load models
if ~exist('cnan','var')
    startcna(1)
end
%% Example 1
% cnap = CNAsbmlModel2MFNetwork('C:\Users\phili\Dokumente\Python\mcs\examples\SmallExample.sbml');
% cnap.reacID = cnap.reacID(:,3:end);
% clear modules
% modules{1}.type = 'lin_constraints';
% modules{1}.sense = 'target';
% modules{1}.V = -ismember(cellstr(cnap.reacID),'R04')';
% modules{1}.v = -1;
% modules{1}.V(2,:) = 2*-ismember(cellstr(cnap.reacID),'R01')';
% modules{1}.v(2,1) = -0.5;
% modules{2}.type = 'lin_constraints';
% modules{2}.sense = 'desired';
% modules{2}.V = -ismember(cellstr(cnap.reacID),'R03')';
% modules{2}.v = -1;
% 
% koCost = [nan nan nan nan 2 3 4 4 4 4]';
% kiCost = [nan 6 nan nan nan nan nan nan nan nan]';
% maxCost = 50;
% maxSolutions = 50;
% options.milp_bigM = 0;
% options.milp_time_limit = inf;
% 
% obj = MCS_cplex(cnap,modules,koCost,kiCost,maxCost,options.milp_bigM);
% [mcs,status] = obj.findMCS(maxSolutions,options.milp_time_limit);

%% Example 2
cnap = CNAsbmlModel2MFNetwork('C:\Users\phili\Dokumente\Python\mcs\examples\weak_coupling.sbml');
clear modules
modules{1}.type = 'bilev_w_constr';
modules{1}.sense = 'desired';
modules{1}.c = zeros(1,cnap.numr);
modules{1}.c(ismember(cellstr(cnap.reacID),'r_BM')) = -1;
modules{1}.V = -ismember(cellstr(cnap.reacID),'r_P')';
modules{1}.v = -1;
modules{1}.V(2,:) = -ismember(cellstr(cnap.reacID),'r_BM')';
modules{1}.v(2,1) = -0.5;

koCost = nan(cnap.numr,1);
koCost(ismember(cellstr(cnap.reacID),{'r2','r3','r4','r5','r6','r7','r8','r9'})) = 1;
kiCost = nan(cnap.numr,1);
kiCost(ismember(cellstr(cnap.reacID),'r1')) = 1;

maxCost = 50;
maxSolutions = 50;
options.milp_bigM = 0;
options.milp_time_limit = inf;

obj = MCS_cplex(cnap,modules,koCost,kiCost,maxCost,options.milp_bigM);
[mcs,status] = obj.findMCS(maxSolutions,options.milp_time_limit);