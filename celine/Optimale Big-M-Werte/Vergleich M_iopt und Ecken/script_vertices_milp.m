%% Script to compre Big-M values derived from vertex enumeration and MILP

clear all
%% 1. select network

% Tiny net
%load('tinynet.mat'); %Target: r_B_ex,
%disp('Tiny net');

% Small Example 3
load('smallex3.mat'); %Target: r_BM,
disp('Tiny net');

%% 3. Build Polyhedron (assume all are knockable)
n = size(A_ineq,1);
m = size(A_ineq,2);
idx_S=(m-n+1:m);

%add farkas-constraint
A_ineq = [A_ineq; A_eq(1,:)];
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
M_ineq   = Max(idx_S(1):idx_S(end));

%% MILP
disp('Computing best bounding M values using MILP.');

rankA=rank(A);
m = size(A,1);
n = size(A,2);
itermax = 1;
T=zeros(1,itermax);

% find random b which gurantees linear independence
A_lin = A';
b_lin = zeros(size(A_lin,1),1);
b_lin(any(A_lin,2)) = 20*(rand(sum(any(A_lin,2)),1) + 1)-30;
% use MILP to ensure supp(y)>=n, minimize the number of y-Variables y \neq 0
% to get a linear combination of b_lin with rows of A
cpx1.f     = [zeros(1,size(A_lin,2)) ones(1,size(A_lin,2))];
cpx1.sense = 'minimize';
cpx1.Aineq = zeros(0,2*size(A_lin,2));
cpx1.bineq = zeros(0,1);
cpx1.Aeq   = [A_lin zeros(size(A_lin,1),size(A_lin,2))];
cpx1.beq   = b_lin;
cpx1.lb    = [-inf(size(A_lin,2),1); zeros(size(A_lin,2),1)];
cpx1.ub    = [ inf(size(A_lin,2),1); ones( size(A_lin,2),1)];
cpx1.ctype = [repmat('C',1,size(A_lin,2)) repmat('B',1,size(A_lin,2))];
cpx1 = Cplex(cpx1); % Create Cplex class

cpx1.addIndicators( (size(A_lin,2)+(1:size(A_lin,2)))',...
    ones(1,size(A_lin,2)),...
    num2cell([eye(size(A_lin,2)) zeros(size(A_lin,2),size(A_lin,2))]',1),...
    repmat('E',size(A_lin,2),1),...
    zeros(size(A_lin,2),1)',...
    '',...
    ones(1,size(A_lin,2)),...
    0.0*ones(size(A_lin,2),1));
cpx1.DisplayFunc = [];
disp('Finding a suitable vector for guaranteeing linear independence.');
cpx1.solve()
while cpx1.Solution.objval < rankA
    disp('Trying again.');
    b_lin(any(A_lin,2)) = 20*(rand(sum(any(A_lin,2)),1) + 1)-30;
    cpx1.Model.lhs(:) = -inf;
    cpx1.Model.rhs = b_lin;
    cpx1.Model.lhs = b_lin;
    cpx1.solve();
end


for k = 1:itermax %loop to compute average time
    
    %construct cplex object
    
    prob.Aineq = [A zeros(m,m+m)]; %xyz
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
    
    M_ineq2 = zeros(size(idx_S,2),1);
    
    %compute M_i^{opt}
    tic
    for i=1:length(idx_S)
        cpx.Model.obj(idx_S(i))=-1;
        cpx.solve;
        M_ineq2(i) = -cpx.Solution.objval;
        cpx.Model.obj(idx_S(i))=0;
    end
    T(k)=toc; %Time for one whole iteration
    
    %% Compare M_i^opt from MILP and Vertex Enumeration
    j   = 0 ;
    eps = 10^-3;
    
    for i=1:length(M_ineq)
        if M_ineq(i) - M_ineq2(i) < eps
            j = j+1; % j increases if both values are equal
        end
    end
    
    if j == length(M_ineq)
        disp('All M_i derived from MILP computation are valid');
    else
        if j==1
            X = sprintf('1 M_i derived from MILP computation is valid');
        else
            X = sprintf('%d M_i derived from MILP computation are valid', j);
        end
        disp(X)
    end
end
disp(['new Big-M Values:'])
M_ineq
time=sum(T)/itermax;
disp(['Average time:' num2str(time)]);