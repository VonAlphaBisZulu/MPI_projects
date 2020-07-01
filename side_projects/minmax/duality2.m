load('duality_model.mat');
cnap = model;
cnap.reacMin(14)=-5;
cnap.macroDefault = [];
model = CNAcna2cobra(cnap);
selectedRxnList = cellstr(cnap.reacID);
options.targetRxn = 'r_P';
options.numDel = 2;
% [optKnockSol, bilevelMILPproblem] = OptKnock(model, selectedRxnList, options);
% [optKnockSol, bilevelMILPproblem] = OptKnock(model, selectedRxnList, options, constrOpt, prevSolutions, verbFlag, solutionFileNameTmp);

%% OptGene

load('duality_model.mat');
cnap = model;
cnap.reacMin(14)=-5;
cnap.macroDefault = [];
model = convertCNAModelToCbModel(cnap);
model.genes = [];


targetRxn = 'r_P';
substrateRxn = 'r_S';
generxnList = cellstr(cnap.reacID);

[x, population, scores, optGeneSol] = optGene(model, targetRxn, substrateRxn, generxnList);

%% RobustKnock

chemicalInd = findStrPos(cnap.reacID,'r_P');
biomassInd = findStrPos(cnap.reacID,'r_BM');
objectiveInd = findStrPos(cnap.reacID,'r_BM');
knockoutNum = 2;
maxW = 789;
constraintsList.reactions = [];

results=robustKnock(model, chemicalInd, biomassInd, objectiveInd, knockoutNum, maxW, 0, 1, constraintsList);