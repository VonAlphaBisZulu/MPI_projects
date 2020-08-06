cnap = CNAloadNetwork({'MTEx2';1},1,1);

idx.S2Ex  = findStrPos(cnap.reacID,'Ex_B');
idx.S2Up   = findStrPos(cnap.reacID,'Ex_BIn');
idx.growth = findStrPos(cnap.reacID,'Ex_BM');
idx.atpm   = findStrPos(cnap.reacID,'ZM');
idx.S1    = findStrPos(cnap.reacID,'Ex_S');
idx.prod   = findStrPos(cnap.reacID,'Ex_P');
idx.rDes   = findStrPos(cnap.reacID,'R13');
idx.rDes2   = findStrPos(cnap.reacID,'R02');
idx.byprod = findStrPos(cnap.reacID,'Ex_E');

% cnap.reacMin(idx.S2Up) = -30;

%% T1 and T2
minYield = 0.3; % change to 0.4 or 0.3
Cp = 4;
Cs = 6;
Cb = 2;
T1 = full(sparse([ 1          1              1              2            3        ], ...
                 [idx.prod    idx.S1        idx.S2Up        idx.S2Ex    idx.atpm  ], ...
                 [Cp          minYield*Cs   minYield*Cb     1           -1        ],3,cnap.numr));
t1 =  [  0  ; 0 ; -1.5];
T2 = full(sparse([ 1          1                             2         ], ...
                 [idx.prod    idx.S1                        idx.atpm  ], ...
                 [Cp          minYield*Cs                   -1        ],2,cnap.numr));
t2 =  [  0   ; -1.5];
%% D
% Desired: Growth 0.1 at ATPM = 1.5
D1 = full(sparse( [1          2 ], ...
                 [idx.growth idx.atpm], ...
                 [-1         -1],2,cnap.numr));
d1 = [-1 ; -1.5];
% D2 = full(sparse( [1          2 ], ...
%                  [idx.rDes  idx.atpm], ...
%                  [1         -1],2,cnap.numr));
% d2 = [-0.2 ; -1.5];
D2 = full(sparse( [1          2 ], ...
                 [idx.rDes2  idx.atpm], ...
                 [-1         -1],2,cnap.numr));
d2 = [-0.2 ; -1.5];

export_reacs = cell2mat(cellfun(@(x) contains(x,'Ex'),cellstr(cnap.reacID),'UniformOutput',false));
notknockable = export_reacs;

koCost = ones(cnap.numr,1);   % all interventions are weighted as 1
koCost(3) = 0.6;
koCost(4) = 1.6;
koCost(5) = 0.9;
koCost(6) = 0.8;
koCost([17 18 19]) = 0;

knockin_ids             = idx.S2Up;
notknockable(idx.S2Up)  = 0;         % The Knock-In-Reaction MUST be set 'knockable'
                                 % in order to work.

max_num_interv = 4;
time_limit = inf;
max_solutions = 1e5;

[is1, status, obj] = ...
    CNAmultiImposedConstrISfinder(cnap,{T1,T2},{t1,t2},{D1 D2},{d1 d2},notknockable,knockin_ids,koCost,max_num_interv,time_limit,max_solutions);

[is2, status, obj] = ...
    CNAmultiImposedConstrISfinder(cnap,{T1,T2},{t1,t2},{D1 D2},{d1 d2},notknockable,knockin_ids,koCost,max_num_interv,time_limit,max_solutions,0);

%% postprocessing
cnapBU = cnap;

%% validation
for i = 1:size(is,2)
    cnap = cnapBU;
    cnap.reacMin(is(:,i)==-1 | isnan(is(:,i))) = 0;
    cnap.reacMax(is(:,i)==-1 | isnan(is(:,i))) = 0;

    if CNAcheckDesired(cnap,T1,t1)
    disp(['Target feasible in mcs ' num2str(i)]); else
    disp(['Target infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckDesired(cnap,T2,t2)
    disp(['Target feasible in mcs ' num2str(i)]); else
    disp(['Target infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckDesired(cnap,D1,d1)
    disp(['Desired feasible in mcs ' num2str(i)]); else
    disp(['Desired infeasible in mcs ' num2str(i)]); 
    end
    if CNAcheckDesired(cnap,D2,d2)
    disp(['Desired feasible in mcs ' num2str(i)]); else
    disp(['Desired infeasible in mcs ' num2str(i)]); 
    end
%   % FVA
%     [lb, ub] = CNAfluxVariability(cnap);
%     disp([cellstr(cnap.reacID) num2cell([lb ub])]);
    
    % special case: For yield minimization, deactivate sink reaction when source reaction is
    % activated
    if ~isnan(is(idx.S2Up,i))
        is(idx.S2Ex,i*~isnan(is(idx.S2Up,i))) = -0.5;
        cnap.reacMin(idx.S2Ex) = 0;
        cnap.reacMax(idx.S2Ex) = 0;
    end
    
    [Ymin,~,~,~] = CNAoptimizeYield(cnap,-iv(cnap.numr,idx.prod)',-iv(cnap.numr,idx.S1)');
    disp(['Ymin (glc only M/M) = ' num2str(Ymin)]); 
    [Ymin,~,~,~] = CNAoptimizeYield(cnap,-Cp*iv(cnap.numr,idx.prod)',-Cs*iv(cnap.numr,idx.S1)'-Cb*iv(cnap.numr,idx.S2Up)');
    disp(['Ymin (total C balance) = ' num2str(Ymin) 10 '______']); 
end