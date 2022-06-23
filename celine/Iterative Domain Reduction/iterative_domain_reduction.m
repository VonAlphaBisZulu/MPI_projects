%% Iterative Domain Reduction

clear all
%% 1. select network
%Tiny net
% load('tinynet.mat'); %Target: r_B_ex, %use data of computeM input after  FVA procedure in computeM-Function
% obj.M   = 1e1;
% obj.M   = 1e2;
% maxSize = 2; 
% disp('Tiny net');

% Small Example 3
load('smallex3.mat'); %Target: r_BM, 
obj.M = 1e3;
obj.M = 1e4;
maxSize = 3; 
disp('Small Example 3');

omega=maxSize; %upper bound of objective function
n = size(A_ineq,1); 
m = size(A_ineq,2);

% A_ineq has already columns for binary variables

for i = 1:n    %add M and set lower and upper bounds for binaries
A_ineq(i,m-n+i) = -obj.M;
lb(m-n+i)=0;
ub(m-n+i)=1;
end

c = [zeros(m-n,1); ones(n,1)]; %add upper objective bound 
b_ineq = [b_ineq; omega];
A_ineq = [A_ineq; c'];

intcon=(m-n+1:m);
f = zeros(1,m);

options = optimoptions('intlinprog','Display','off','MaxNodes', 1e4); %nodelimit 10k

%split matrix
A_ineq_pos=zeros(size(A_ineq,1),size(A_ineq,2));
A_ineq_neg=zeros(size(A_ineq,1),size(A_ineq,2));
             for j=1:size(A_ineq,2)
                 for i=1:size(A_ineq,1)
                     if A_ineq(i,j)>=0 
                         A_ineq_pos(i,j)=A_ineq(i,j); %only positive entries
                     else
                         A_ineq_neg(i,j)=A_ineq(i,j); %only negative entries
                     end
                 end
             end 

count = 0;  % count for break up criteria
iter  = 0;  % number of iterations
tol   = 0.1;% tolerance
while count < n
 count = 0; 
 iter = iter +1;
for j=1:(m-n)

f(j) = -1;     %maximize continuous variables
[v,fval] = intlinprog(f,intcon,A_ineq,b_ineq,A_eq,b_eq,lb,ub,[], options); 
ub(j) = -fval; %save as new upper bound
f(j) = 0;

f(j) = 1;      %minimize continuous variables
[v,fval] = intlinprog(f,intcon,A_ineq,b_ineq,A_eq,b_eq,lb,ub,[], options);
f(j) = 0;
lb(j) = fval;  %save as new lower bound

end
%compute new M-Values
for i=1:n
M(i)=dot(ub(1:(m-n)), A_ineq_pos(i,1:(m-n)))+dot(lb(1:(m-n)), A_ineq_neg(i,1:(m-n)));
if abs(A_ineq(i,m-n+i)+M(i)) < tol
    count=count+1;
end
A_ineq(i,m-n+i)=-M(i);
end

end

%display output
disp(['M-Values:']);
M'
disp(['new lower bounds:']);
lb(1:m-n)
disp(['new upper bounds:']);
ub(1:m-n)
disp(['Iterations:' num2str(iter-1)]);

