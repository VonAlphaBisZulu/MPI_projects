function [M_ineq, b_eq0 ] = computeM(obj,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_ineq) %ECC2Comp -10^0
%select b_F
b_eq0 = -10^(0);


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

%%%%%%%%%%%%%% Approach Scaling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(A_ineq,1);
m = size(A_ineq,2);
lb(isinf(lb))=-100;
ub(isinf(ub))=100;

%split Matrix
A_ineq_pos=zeros(size(A_ineq,1),size(A_ineq,2));
A_ineq_neg=zeros(size(A_ineq,1),size(A_ineq,2));
for j=1:size(A_ineq,2)
    for i=1:size(A_ineq,1)
        if A_ineq(i,j)>=0
            A_ineq_pos(i,j)=A_ineq(i,j); %only positiv entries
        else
            A_ineq_neg(i,j)=A_ineq(i,j); %only negativ entries
        end
    end
end
%compute M
for i=1:n
    M_ineq(i) = dot(ub(1:m), A_ineq_pos(i,1:m))+dot(lb(1:(m)), A_ineq_neg(i,1:(m)));
end