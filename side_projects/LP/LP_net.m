startcna(1);
model = 'SmallExample2';
cnap = CNAloadNetwork(find(strcmp(cellstr({cnan.net.path}'),model)),1,1);

%% Modify model and constrains
% Limit uptake of A and B, min Reaction Rate of E export, "negative" export
% rates
cnap.reacMin(1) = -1; % maximum uptake of A is 1
cnap.reacMax(1) =  0;
cnap.reacMin(8) = -1; % maximum uptake/excretion of B is 1
cnap.reacMax(8) =  1;
cnap.reacMin(10) = 0.0; % minimum rate for export of E
cnap.stoichMat(1,1) = -1; % Transport over system bounds is defined negatively
cnap.stoichMat(2,8) = -1; % Transport over system bounds is defined negatively

cnap.stoichMat(:,1) = cnap.stoichMat(:,1)*1.3; % Transport over system bounds is defined negatively
cnap.stoichMat(:,3) = cnap.stoichMat(:,3)*0.33;
cnap.stoichMat(:,4) = cnap.stoichMat(:,4)*1.2;
cnap.stoichMat(:,8) = cnap.stoichMat(:,8)*2.7;
cnap.stoichMat(:,9) = cnap.stoichMat(:,9)*0.9;
cnap.stoichMat(:,10) = cnap.stoichMat(:,10)*1.2;

% cnap.objFunc(5) = 1; % optimize for P export

%% Map reversible reactions as two IRREV and update stoichiometric Matrix
cnapIrr = cnap;

cnapIrr.stoichMat(:,end+1:2*end) = -cnap.stoichMat; % add negative stoichiometric
                                                    % Matrix for reverse sense of
                                                    % reactions to the end
                                                    % of the st. Matrix
cnapIrr.reacMax(end+1:2*end) = -cnap.reacMin;
cnapIrr.reacMin(end+1:2*end) = -cnap.reacMax;       % update boundaries
cnapIrr.reacMin(cnapIrr.reacMin<0) = 0;
cnapIrr.reacMax(cnapIrr.reacMax<0) = 0;
% add identifiers
cnapIrr.reacID = char([cellstr(cnap.reacID(:,:)); strcat(cellstr(cnap.reacID(:,:)),'_rev')]);

% make high upper bound -> inf
cnapIrr.reacMax(cnapIrr.reacMax>=100) = Inf;

numReacs = size(cnapIrr.stoichMat,2);
numMets = size(cnapIrr.stoichMat,1);
varIndices = repmat({'rp'},numReacs/2,1); % reactions in positive direction
varIndices = [varIndices; repmat({'ri'},numReacs/2,1)]; % reactions in inverse direction

% Check if stoichiometric network is feasible
fba = Cplex();
fbaOpts = cplexoptimset;
fba.Model.sense = 'maximize';
fba.Model.A = cnapIrr.stoichMat;
fba.Model.obj = zeros(size(cnapIrr.stoichMat,2),1);
fba.Model.lb = cnapIrr.reacMin;
fba.Model.ub = cnapIrr.reacMax;
fba.Model.lhs = zeros(numMets, 1);
fba.Model.rhs = zeros(numMets, 1);
refConf = fba.solve();

%% ======== EM calculation =========

% Create A from IRREVersible stoichio. Matrix
A = cnapIrr.stoichMat;
% define and add indicator values
% add rows, where reaction and reverse are "added" to verify that there is
% only one active
A(end+1:end+0.5*numReacs,end+1:(end+numReacs)) = [eye(numReacs/2) eye(numReacs/2)]; %
varIndices = [varIndices; repmat({'idc'},numReacs,1)]; % indicator variables
ctype = [repmat('C',1,numReacs) repmat('B',1,numReacs)];

clear indicators

%% Defining indicator Values
for i = 1:numReacs
    indicators(i).variable      = numReacs+i;
    indicators(i).complemented	= 1  ;
    a   = zeros(2*numReacs, 1);
    a(i)= 1.0;
    indicators(i).a             = a  ;
    indicators(i).sense         = 'E';
    indicators(i).rhs           = 0.0 ;
    indicators(i).name          = cellstr(['ind_' cnapIrr.reacID(i,:)]);
end

% Make row that says: sum of all indicator values (in other words: active reactions)

A(end+1,:) = strcmp(varIndices,'idc')';
A(end+1,:) = [ones(numReacs, 1); zeros(numReacs, 1)]; % exclude trivial solution
        
cpx = Cplex();
cpxOpts = cplexoptimset;
cpx.Model.sense = 'minimize';
cpx.Model.A = A;
cpx.Model.obj = double(strcmp(varIndices,'idc'));
for i =1:numReacs % add structure doesn't work for some reason
    cpx.addIndicators(indicators(i).variable,...
                  indicators(i).complemented,...
                  indicators(i).a,...
                  indicators(i).sense,...
                  indicators(i).rhs,...
                  indicators(i).name);
end
cpx.Model.lb = [cnapIrr.reacMin; zeros(numReacs,1)];
cpx.Model.ub = [cnapIrr.reacMax; ones(numReacs,1)];
cpx.Model.lhs = [zeros(numMets, 1); zeros(0.5*numReacs, 1); 6;1];
cpx.Model.rhs = [zeros(numMets, 1); ones(0.5*numReacs, 1); 0.5*numReacs;Inf];
cpx.Model.ctype = ctype;
sol = [];
sol = cpx.populate();

if sol.status == 103
    cpx.refineConflict();
end
try
    clear solutions
    for i=1:length(cpx.MipStart)
        solutions(i).reacs = [cellstr(cnapIrr.reacID(cpx.MipStart(i).x(1:20)~=0,:)),...
        num2cell(cpx.MipStart(i).x(cpx.MipStart(i).x(1:20)~=0))];
        ind = find(cpx.MipStart(i).x(1:40)~=0);
        solutions(i).ind = [cellstr(cnapIrr.reacID(cpx.MipStart(i).x(21:40)~=0,:)),...
        num2cell(cpx.MipStart(i).x(ind(ind>=21)))];
    end
    for i=1:length(cpx.MipStart)
        if length(solutions(i).ind)==round(sum(sol.x(21:40)))
            if length(solutions(i).ind)==length(solutions(i).reacs)
                disp(solutions(i).reacs)
            end
        end
    end
catch
end

%% FBA
%{
cpx = Cplex();
cpxOpts = cplexoptimset;
cpx.Model.sense = 'maximize';
cpx.Model.obj = cnapIrr.objFunc;
cpx.Model.lb = cnapIrr.reacMin;
cpx.Model.ub = cnapIrr.reacMax;
cpx.Model.A = cnapIrr.stoichMat;
cpx.Model.lhs = zeros(size(cnapIrr.stoichMat, 1), 1);
cpx.Model.rhs = zeros(size(cnapIrr.stoichMat, 1), 1);
sol = cpx.solve();

% Map solution to original network
solution = repmat({'' NaN},[size(cnap.stoichMat,2),1]);
for i = 1:size(cnap.stoichMat,2)
    solution{i,1} = cnap.reacID(i,:);
    solution{i,2} = sol.x(i);
    if any(mapIrr(:,1)==i)
        if sol.x(mapIrr(mapIrr(:,1)==i,2)) > 0
            solution{i,2} = -sol.x(mapIrr(mapIrr(:,1)==i,2));
        end
    end
end
solution
%}
