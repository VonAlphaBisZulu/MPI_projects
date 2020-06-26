cnap = robustknock;
% translating
A  = cnap.stoichMat;
b  = zeros(cnap.nums,1);
lb = cnap.reacMin;
ub = cnap.reacMax;
c  = cnap.objFunc;
n  = cnap.numr;
% splitting reversible reactions
A  = [A, -A(:,lb < 0)];
ub = [ub; -lb(lb < 0)];
c  = [c; -c(lb < 0)];
% transform to min(c'*x) s.t. A*x <= 0
A = [A;-A;eye(n)]; % to mime A*x = 0 and upper bound
b = [b;-b;ub];

%% solve primal
[o,q,p] = cplexlp(c,A,b,[],[],zeros(length(c)),inf(length(c)));

%% solve dual
[o2,q2,p2] = cplexlp(b,-A',c,[],[],zeros(length(b)),inf(length(b)));