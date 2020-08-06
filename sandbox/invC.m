if ~exist('cnan','var')
    startcna(1)
end
cnap = CNAloadNetwork({'mtest';1},1,1);
cnapBU = cnap;

% r0  = find(strcmp(cellstr(cnap.reacID),'rz'));
% cnap.reacMin(r0) = 0;
% cnap.reacMax(r0) = 0;
% r0  = find(strcmp(cellstr(cnap.reacID),'r1'));
% cnap.reacMin(r0) = 0;
% cnap.reacMax(r0) = 0;


% indices
p_ex = find(strcmp(cellstr(cnap.reacID),'rp_ex'));
t_ex = find(strcmp(cellstr(cnap.reacID),'rt_ex'));
t_up = find(strcmp(cellstr(cnap.reacID),'rt_up'));
s_up = find(strcmp(cellstr(cnap.reacID),'rs_up'));
bm   = find(strcmp(cellstr(cnap.reacID),'rbm'));
z    = find(strcmp(cellstr(cnap.reacID),'rz'));

ft = -1;
fs = -1;
Y  = 0.4;

%% Target and Desired Polytopes
% Target minimum Yield of 
T1 = full(sparse(  [1     1     2     3], ...
                   [p_ex  s_up  t_up  bm  ], ...
                   [1     fs*Y  1     -1  ],3,cnap.numr));
t1 =  [  0 ; 0 ; -1];
T2 = full(sparse(  [1     1      1      2      3  ], ...
                   [p_ex  s_up   t_up   t_ex   bm  ], ...
                   [1     fs*Y   ft*Y   1      -1  ],3,cnap.numr));
t2 =  [  0 ; 0 ; -1];
% T2 = T1;
% t2 = t1;
% T1 = T2;
% t1 = t2;
% [a,b,c,e] = CNAoptimizeFlux(cnap,[],[],2, 0, 0, T2, t2);
% [a,b,c,e] = CNAoptimizeYield(cnap,[],[],2);

% Desired: Growth 0.1 at maintanance 0.1
D1 = full(sparse( [1      2 ], ...
                  [bm     z ], ...
                  [-1    -1 ],2,cnap.numr));
d1 = [-1 ; -1];
D2 = full(sparse( [1      2 ], ...
                  [p_ex   z ], ...
                  [-1    -1 ],2,cnap.numr));
d2 = [-4 ; -1];

% Knockout Costs
cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',repmat({''},1,cnap.numr));
[cnap,~,genes,gr_rules] = CNAgenerateGERassociation(cnap);

rkoCost          = ones(cnap.numr,1);
rkoCost([p_ex t_ex t_up s_up bm z]) = nan;
% rkoCost          = nan(cnap.numr,1);
rki  = findStrPos(cnap.reacID,'r[5,6]|rt_up','regex');
rkiCost      = nan(cnap.numr,1);
rkiCost(rki) = 1;

gkoCost                      = nan(length(genes),1);
% gkoCost                      = ones(length(genes),1);
% gkoCost(strcmp(genes,'sp'))  = nan;
gkiCost                      = nan(length(genes),1);
% gkiCost(findStrPos(genes,{'a1' 'a2' 'a3'})) = 0.15; % KI of EX etoh

%% computation

% MCS search parameters
max_num_interv  = inf;
max_solutions = inf;
time_limit = inf;
default_flux_limit = 1000;

[gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx] = CNAgeneMultiImposedMCSFinder(cnap, {T1,T2} , {t1,t2} , {D1,D2} , {d1,d2} ,...
                                                rkoCost,rkiCost, ...
                                                max_num_interv,time_limit,max_solutions,...
                                                1,0, ... use_compression,use_bigM,
                                                1,gkoCost,gkiCost); % enum_method, gkoCost, gkiCost
