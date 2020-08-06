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

numReacs = size(cnap.stoichMat,2);
numMets = size(cnap.stoichMat,1);
varIndices = repmat({'r'},numReacs,1); % reactions in positive direction

cnap.reacMin(cnap.reacMin<=-100) = -Inf;
cnap.reacMax(cnap.reacMin>= 100) =  Inf;
obj=zeros(size(cnap.stoichMat,2),1);
obj(5)=1;

% Check if stoichiometric network is feasible
fba = Cplex();
fbaOpts = cplexoptimset;
fba.Model.sense = 'maximize';
fba.Model.A = cnap.stoichMat;
fba.Model.obj = obj;
fba.Model.lb = cnap.reacMin;
fba.Model.ub = cnap.reacMax;
fba.Model.lhs = zeros(numMets, 1);
fba.Model.rhs = zeros(numMets, 1);
refConf = fba.solve();

%% ======== EM calculation =========

% Create A from IRREVersible stoichio. Matrix
A = cnap.stoichMat;
% define and add indicator values
varIndices = [varIndices; repmat({'idc+'},numReacs,1); repmat({'idc-'},numReacs,1)];
ctype = [repmat('C',1,numReacs) repmat('B',1,2*numReacs)];
A(end+1:end+numReacs,1:end+2*numReacs) = [ zeros(numReacs) eye(numReacs) eye(numReacs)  ];

clear indicators

%% Defining indicator Values
for i = 1:numReacs
    indicators(i).variable      = numReacs+i;
    indicators(i).complemented	= 0  ;
    indicators(i).a             = iv(3*numReacs,i);
    indicators(i).sense         = 'E';
    indicators(i).rhs           = -1e-10 ;
    indicators(i).name          = cellstr(['ind+_' cnap.reacID(i,:)]);
end
for i = 1:numReacs
    indicators(i+numReacs).variable      = numReacs+i;
    indicators(i+numReacs).complemented	 = 1  ;
    indicators(i+numReacs).a             = iv(3*numReacs,i);
    indicators(i+numReacs).sense         = 'L';
    indicators(i+numReacs).rhs           = 1e-10;
    indicators(i+numReacs).name          = cellstr(['ind-_' cnap.reacID(i,:)]);
end

% Make row that says: sum of all indicator values (in other words: active reactions)

A(end+1,:) = strcmp(varIndices,'idc+')' | strcmp(varIndices,'idc-')';
%A(end+1,:) = [ones(numReacs, 1); zeros(numReacs, 1)]; % exclude trivial solution
        
cpx = Cplex();
cpxOpts = cplexoptimset;
cpx.Model.sense = 'minimize';
cpx.Model.A = A;
cpx.Model.obj = double(strcmp(varIndices,'idc+')|strcmp(varIndices,'idc-'));
for i =1:2*numReacs % add structure doesn't work for some reason
    cpx.addIndicators(indicators(i).variable,...
                  indicators(i).complemented,...
                  indicators(i).a,...
                  indicators(i).sense,...
                  indicators(i).rhs,...
                  indicators(i).name);
end

cpx.Model.lb = [cnap.reacMin; zeros(2*numReacs,1)];
cpx.Model.ub = [cnap.reacMax; ones(numReacs,1) ; zeros(numReacs,1)];
cpx.Model.lhs = [zeros(numMets, 1); zeros(numReacs,1); 1];
cpx.Model.rhs = [zeros(numMets, 1); ones(numReacs,1); Inf];
cpx.Model.ctype = ctype;
sol = [];
sol = cpx.populate();

if sol.status == 103
    cpx.refineConflict();
end
try
    clear solutions
    for i=1:length(cpx.MipStart)
        solutions(i).reacs = [cellstr(cnapIrr.reacID(cpx.MipStart(i).x(1:numReacs)~=0,:)),...
        num2cell(cpx.MipStart(i).x(cpx.MipStart(i).x(1:numReacs)~=0))];
        ind = find(cpx.MipStart(i).x(1:2*numReacs)~=0);
        solutions(i).ind = [cellstr(cnapIrr.reacID(cpx.MipStart(i).x(numReacs+1:2*numReacs)~=0,:)),...
        num2cell(cpx.MipStart(i).x(ind(ind>=numReacs+1)))];
    end
    for i=1:length(cpx.MipStart)
        if length(solutions(i).ind)==round(sum(sol.x(numReacs+1:2*numReacs)))
            if length(solutions(i).ind)==length(solutions(i).reacs)
                disp(solutions(i).reacs)
            end
        end
    end
catch
end