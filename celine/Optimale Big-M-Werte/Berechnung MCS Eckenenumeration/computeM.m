function [M_ineq, b_eq0] = computeM(obj,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_ineq) %function for vertex approach
b_eq0=b_eq;

% if any variable is fixed, add it to b
fixed_vars = (lb == ub);
A_ineq = A_ineq(:,~fixed_vars);
A_eq   = A_eq(:,~fixed_vars);
b_ineq = b_ineq - sum(A_ineq(:,fixed_vars).*ub(fixed_vars)',2);
b_eq   = b_eq - sum(A_eq(:,fixed_vars).*ub(fixed_vars)',2);
lb = lb(~fixed_vars);
ub = ub(~fixed_vars);

%% This function is used to compute bigM for each constraint/bound
ko_ineq = any(z_ineq,1)';
M_ineq = nan(size(A_ineq,1),1);
%% Find boundaries on constraints via direct optimization (FVA-like)

for i = find(ko_ineq)'
    lp.c = -A_ineq(i,:);
    lp.A_ineq = A_ineq(~ko_ineq,:);
    lp.b_ineq = b_ineq(~ko_ineq,:);
    lp.A_eq = A_eq;
    lp.b_eq = b_eq;
    lp.lb = lb;
    lp.ub = ub;
    x = obj.solveLP(lp);
    if ~isempty(x)
        disp(['Found bound ' num2str(-lp.c*x)]);
        M_ineq(i) = -lp.c*x+100; % maybe + epsilon?
    else
        M_ineq(i) = obj.M;
    end
    
end

foundM = ~isnan(M_ineq) & ~isinf(M_ineq);
ko_ineq(foundM) = false;
A_ineq = [A_ineq   full(sparse(find(foundM),1:sum(foundM),-1,size(A_ineq,1),sum(foundM)))];
A_eq   = [A_eq    zeros(size(A_eq,1),sum(foundM))];
lb     = [lb; zeros(sum(foundM),1)];
ub     = [ub; M_ineq(foundM)];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Vertex Computation %%%%%%%%%%%%%%%%%%%%%%%%
n = size(A_ineq,1);
m = size(A_ineq,2);
idx_S=(m-n+1:m);
A_ineq = [A_ineq; A_eq(1,:)]; %add farkas, computation only for Target Region!
b_ineq = [b_ineq; b_eq(1,:)];
lb_bounded=~isinf(lb);
% add lower bounds
LB = full(sparse(1:sum(lb_bounded),find(lb_bounded),-1,sum(lb_bounded),size(A_ineq,2)));
A=[A_ineq; LB];
b=[b_ineq; -lb(lb_bounded)];

%% 4. Compute Vertices via MPT3-Tool

P1=Polyhedron('A', A, 'b', b);
P1.computeVRep
P1.V;
%select M_i opt
Max = max(P1.V,[], 1)';
M_ineq   =Max(idx_S(1):idx_S(end));
disp(['Big-M Values which were found by using vertex enumeration :'])
M_ineq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%