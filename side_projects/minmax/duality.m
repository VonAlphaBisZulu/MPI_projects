load('duality_model.mat');
cnap = model;
cnap.reacMin(ismember(cellstr(cnap.reacID),'r_6')) = 0;
cnap.reacMax(ismember(cellstr(cnap.reacID),'r_6')) = 0;
% translating
A  = cnap.stoichMat;
b  = zeros(cnap.nums,1);
lb = cnap.reacMin;
ub = cnap.reacMax;
c  = cnap.objFunc;
n  = size(cnap.stoichMat,2);

%% 0. Test original model with cplexLP

[~,~,~,opt_ref] = CNAoptimizeFlux(cnap);

% opt_p0 is the optimal value in the primal problem
[~,opt_p0,~] = cplexlp(c,[],[],A,b,lb,ub);

if isempty(opt_p0) || abs(opt_ref-opt_p0)>1e-7
    error('LP didn''t return correct solution');
end

%% 1. Split all reactions, optimize with cplexLP and verify if identical solution is found
A_p1  = [A  , -A ];
ub1 = max(zeros(2*n,1),[ub;-lb]);
lb1 = max(zeros(2*n,1),[lb;-ub]);
c_p1 = [c;-c];

[~,opt_p1,~] = cplexlp(c_p1,[],[],A_p1,b,lb1,ub1);

if isempty(opt_p1) || abs(opt_ref-opt_p1)>1e-7
    error('LP didn''t return correct solution');
end

%% 2. Reshape this problem into A*x<=0, x>=0 and min(c'*x)
A_p2 = [A_p1 ;-A_p1;eye(2*n);-eye(2*n)]; % to mime A*x = 0 and upper bound
b_p2 = [ b ; b ;  ub1   ;  -lb1   ];
c_p2 = c_p1;

[~,opt_p2,~] = cplexlp(c_p2,A_p2,b_p2,[],[],zeros(length(c_p2),1),inf(length(c_p2),1));

if isempty(opt_p2) || abs(opt_ref-opt_p2)>1e-7
    error('LP didn''t return correct solution');
end

%% 3. Solve dual problem
%
% Primal: min(c'x),  A *x<=b, x>=0
% Dual:   max(b'y),  A'*y>=c, y>=0
%
%  or in standard form for minimization:
%
% Dual:  -min(b'y), -A'*y<=c, y>=0
%
%
%
A_d1 = -A_p2';
b_d1 = c_p2;
c_d1 = b_p2;
% 
[~,q_d1,~] = cplexlp(c_d1,A_d1,b_d1,[],[],zeros(length(c_d1),1),inf(length(c_d1),1));

if isempty(q_d1) || abs(opt_ref+q_d1)>1e-7
    error('LP didn''t return correct solution');
end

%% 4. Connect to primal problem
A_ineq_m1 = [zeros(size(A_d1,1),size(A,2)) , A_d1];
b_ineq_m1 = b_d1;
A_eq_m = [A , zeros(size(A,1),size(A_d1,2)) ];
b_eq_m = b;
% adding a constraint to connect both biomass rates
A_eq_m = [A_eq_m ; c' c_d1'];
b_eq_m = [b_eq_m; 0];

lb_m = [lb; zeros(length(c_d1),1)];
ub_m = [ub; inf(length(c_d1),1)];

c_prod = zeros(1,size(A_ineq_m1,2));
c_prod(2) = -1;
[~,q_d2,~] = cplexlp(c_prod,A_ineq_m1,b_ineq_m1,A_eq_m,b_eq_m,lb_m,ub_m);