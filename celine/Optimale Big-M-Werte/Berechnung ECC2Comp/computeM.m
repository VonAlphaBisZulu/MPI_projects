function [M_ineq, b_eq0] = computeM(obj,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_ineq) %function for MILP-M_i^{opt} approach
Offset = 0; % choose one or zero do decide if an offset should be used
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
        M_ineq(i) = -lp.c*x+100;
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

%%%%%%%%%%%% Using MILP to Compute M_I opt %%%%%%%%%%%%%%%%%%%%%%%%

% Build Polyhedron
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

disp('Computing best bounding M values using MILP.');

rankA = rank(A);
m = size(A,1);
n = size(A,2);
itermax = 1;
T=zeros(1,itermax);

% find random b which gurantees linear independence
A_lin = A';
b_lin = zeros(size(A_lin,1),1);
b_lin(any(A_lin,2)) = 20*(rand(sum(any(A_lin,2)),1) + 1)-30;
% use MILP to ensure supp(y)>=n

% while using ECC2comp: it takes to long, in general a random vector does
% the job

% cpx1.f     = [zeros(1,size(A_lin,2)) ones(1,size(A_lin,2))];
%                 cpx1.sense = 'minimize';
%                 cpx1.Aineq = zeros(0,2*size(A_lin,2));
%                 cpx1.bineq = zeros(0,1);
%                 cpx1.Aeq   = [A_lin zeros(size(A_lin,1),size(A_lin,2))];
%                 cpx1.beq   = b_lin;
%                 cpx1.lb    = [-inf(size(A_lin,2),1); zeros(size(A_lin,2),1)];
%                 cpx1.ub    = [ inf(size(A_lin,2),1); ones( size(A_lin,2),1)];
%                 cpx1.ctype = [repmat('C',1,size(A_lin,2)) repmat('B',1,size(A_lin,2))];
%                 cpx1 = Cplex(cpx1); % Create Cplex class
%
%                 cpx1.addIndicators( (size(A_lin,2)+(1:size(A_lin,2)))',...
%                                     ones(1,size(A_lin,2)),...
%                                     num2cell([eye(size(A_lin,2)) zeros(size(A_lin,2),size(A_lin,2))]',1),...
%                                     repmat('E',size(A_lin,2),1),...
%                                     zeros(size(A_lin,2),1)',...
%                                     '',...
%                                     ones(1,size(A_lin,2)),...
%                                     0.0*ones(size(A_lin,2),1));
%                 cpx1.DisplayFunc = [];
%                 disp('Finding a suitable vector for guaranteeing linear independence.');
%                 cpx1.solve()
%                 while cpx1.Solution.objval < rankA
%                     disp('Trying again.');
%                     b_lin(any(A_lin,2)) = 20*(rand(sum(any(A_lin,2)),1) + 1)-30;
%                     cpx1.Model.lhs(:) = -inf;
%                     cpx1.Model.rhs = b_lin;
%                     cpx1.Model.lhs = b_lin;
%                     cpx1.solve();
%                 end


for k = 1:itermax %loop to compute average time
    
    %construct cplex object
    
    prob.Aineq = [A zeros(m,m+m)];
    prob.bineq = b;
    prob.Aeq = [zeros(n,n) A_lin zeros(n,m)];
    prob.beq = b_lin;
    ctype=char(ones([1 (size(prob.Aineq,2))])*('C'));
    ctype(n+m+1:n+m+m)='B';
    prob.ctype = ctype;
    prob.f = zeros(size(prob.Aineq,2),1);
    cpx=Cplex(prob);
    cpx.DisplayFunc = [];
    
    % add indicators part 1 z==1 => a'x=b holds
    indicator.sense = char(ones([1 m])*('E')); %determine the type of the constraints, here: Equality
    indicator.variable = (n+m+1:n+m+m);
    indicator.complemented=zeros(1, size(indicator.variable,2)); %zero: z=1 =>ax=b holds!
    a=cell(1, size(indicator.variable,2)); %coefficients for indicator constraints
    for i=1:size(indicator.variable,2)
        a{i}=prob.Aineq(i,:)';
    end
    indicator.a=a;
    indicator.rhs=b;
    
    cpx.addIndicators(indicator.variable, indicator.complemented, indicator.a, indicator.sense, indicator.rhs);
    
    % add indicators part 2 z==0 => y = 0 holds
    indicator2.sense = char(ones([1 m])*('E'));
    indicator2.variable = (n+m+1:n+m+m);
    indicator2.complemented=ones(1, size(indicator2.variable,2));
    a2=cell(1, size(indicator2.variable,2));
    Ainputind2 = [zeros(m,n) eye(m) zeros(m,m)];
    
    for i=1:size(indicator.variable,2)
        a2{i}=Ainputind2(i,:)';
    end
    indicator2.a=a2;
    indicator2.rhs=zeros(m,1);
    
    cpx.addIndicators(indicator2.variable, indicator2.complemented, indicator2.a, indicator2.sense, indicator2.rhs);
    cpx.Param.mip.limits.nodes.Cur  = 5*10^4; %nodelimit
    M_ineq = zeros(size(idx_S,2),1);
    
    %compute M_i^{opt}
    
    tic
    for i=1:length(idx_S)
        cpx.Model.obj(idx_S(i))=-1;
        cpx.solve;
        M_ineq(i) = -cpx.Solution.objval;
        cpx.Model.obj(idx_S(i))=0;
    end
    T(k)=toc; %Time for one whole iteration
end
disp(['new Big-M Values:'])
M_ineq
time=sum(T)/itermax;
disp(['Time:' num2str(time)]);

if Offset == 1 % Add Offset
    off= 0.69*max(M_ineq);
    M_ineq = M_ineq + off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%