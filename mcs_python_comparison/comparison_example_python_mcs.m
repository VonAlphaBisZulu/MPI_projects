%% Initialize and load models
if ~exist('cnan','var')
    startcna(1)
end
cnap = CNAsbmlModel2MFNetwork('C:\Users\phili\Dokumente\Python\mcs\examples\SmallExample.sbml');
cnap.reacID = cnap.reacID(:,3:end);
modules{1}.type = 'lin_constraints';
modules{1}.sense = 'target';
modules{1}.V = -ismember(cellstr(cnap.reacID),'R4')';
modules{1}.v = -1;
modules{2}.type = 'lin_constraints';
modules{2}.sense = 'desired';
modules{2}.V = -ismember(cellstr(cnap.reacID),'R3')';
modules{2}.v = -1;

koCost = [nan nan nan nan 2 3 4 4 4 4]';
kiCost = [nan 6 nan nan nan nan nan nan nan nan]';
maxCost = 50;
maxSolutions = 50;
options.milp_bigM = 0;
options.milp_time_limit = inf;

obj = MCS_cplex(cnap,modules,koCost,kiCost,maxCost,options.milp_bigM);
[mcs,status] = obj.findMCS(maxSolutions,options.milp_time_limit);