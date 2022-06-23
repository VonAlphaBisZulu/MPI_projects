%% Script AMIP, ALP
clear all
%% 1. select network
%Tiny net
% load('tinynet.mat'); %Target: r_B_ex, %use data of computeM input
% obj.M   = 1e2;
% maxSize = 2;
% disp('Tiny net');

% Small Example 3
load('smallex3.mat'); %Target: r_BM,
obj.M = 1e3;
maxSize = 3;
disp('Small Example 3');


%% methode 1: ALP and AMIP
omega=maxSize;                %value for upper bound of objective function

% construct considerd system: v_i <= M_i*z_i

n = size(A_ineq,1);
m = size(A_ineq,2);
z = [zeros(n,m-n) eye(n) -obj.M*eye(n)];
A_ineq = [A_ineq zeros(n,n)]; %introduce binary variables, A_ineq has already columns for v variables
A_ineq = [A_ineq; z];         %add constraints: v_i <= M_i*z_i
b_ineq = [b_ineq; zeros(n,1)];
A_eq = [A_eq zeros(size(A_eq,1),n)];
lb = [lb; zeros(n,1)];        %add bounds for binaries
ub(m-n+1:m) = inf;
ub = [ub; ones(n,1)];
c = [zeros(m,1); ones(n,1)];

%add constraint: upper bound to objective function

A_ineq = [A_ineq; c']; %add upper bound omega to objective function
A_ineq_store = A_ineq; %save A_ineq
b_ineq = [b_ineq; omega];
f = zeros(m+n,1);
intcon = (m+1: m+n);
options = optimoptions('intlinprog','Display','off');
j = 0;

%% AMIP without timelimit because the systems are small enough
disp('AMIP:')
tol = 0.1;
for k = 1:n %for each M_i one iteration
    f(m-n+k) = -1; %maximize v-Variable
    
    [~,fval] = intlinprog(f,intcon,A_ineq,b_ineq,A_eq,b_eq,lb,ub,[], options);
    if abs(A_ineq(n+k,m+k) - fval) <=tol  %check for changement
        j=j+1;
    end
    f(m-n+k) = 0;
    A_ineq(n+k,m+k) = fval; % put new value in matrix
    
end
if j==n
    disp('No value has changed');
    disp(['Value for M: ' num2str(-fval)]);
else
    disp('Some values have changed');
    for k = 1:n
        disp(['Value for' num2str(k) ': ' num2str(-A_ineq(n+k,m+k))]);
    end
end

%% ALP - linear relaxtion  of AMILP
A_ineq  = A_ineq_store; %reload A_ineq
s       = 1;            %count for outer iteration
eps     = 1e-3;         %tolerance
smax    = 10;           %maximum number of iterations
delta   = inf*ones(smax,1); %initial value for break up criteria
j       = 0;            %changement coutner
options = optimoptions('linprog','Display','off');

disp('ALP:')
while delta(s)> eps && s<smax
    delta_0 = 0;
    for k=1:n
        f(m-n+k) = -1; %maximize v-Variable
        [~,fval] = linprog(f,A_ineq,b_ineq,A_eq,b_eq,lb,ub,[],options);
        delta_0  = delta_0 + A_ineq(n+k,m+k) - fval;
        f(m-n+k) = 0;
        if abs(A_ineq(n+k,m+k) - fval) <=tol
            j=j+1;
        end
        A_ineq(n+k,m+k) = fval;
    end
    s = s+1;
    delta(s) = abs(delta_0);
end

if j==n
    disp('No value has changed');
    disp(['Value for M: ' num2str(-fval)]);
    disp(['Outer Iterations: ' num2str(s-1)]);
else
    disp('Some values have changed');
    disp(['Outer Iterations: ' num2str(s-1)]);
    for k = 1:n
        disp(['Value for' num2str(k) ': ' num2str(-A_ineq(n+k,m+k))]);
    end
end