function [feasible, r] = CNAcheckRegionFeasibility(cnap,D,d,solver)
if nargout >= 2 % do parsimonious FBA if two outputs are requested
    parsimonious = 1;
else
    parsimonious = 0;
end
if nargin < 4
    solver = 0;
end
if nargin == 1
    D = double.empty(0,cnap.numr);
    d = double.empty(0,1);
end
r   = zeros(cnap.numr,1);
N   = initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
ub  = cnap.reacMax;
lb  = cnap.reacMin;
b   = zeros(cnap.numis, 1);
obj = zeros(cnap.numr, 1);
if parsimonious == 1  % Max rates from D vector (workaround for ratios - maximize rates)
    for i = 1:length(d)
        if d(i) ~= 0
            obj(D(i,:)~=0) = -D(i,D(i,:)~=0)*sign(d(i));
        else
            obj(D(i,:)~=0) = -D(i,D(i,:)>0);
            break;
        end
    end
end

switch solver
    case 0 % GLPK       % exitflag ==  2: optimal  % exitflag ==  5: feasible
        ctype = [repmat('S',1,cnap.numis), repmat('U',1,length(d))];
        [x, ~,exitflag] = glpk(obj,[N; D],[b; d],lb,ub,ctype);
    case 1 % linprog    % exitflag ==  1: optimal, feasible
	if ~verLessThan('matlab', '8.1') % use different functions to suppress linprog verbose depending on MATLAB version
    		linprog_options = optimoptions('linprog','Display','off');
	else
    		linprog_options = optimset('Display','off');
	end
        [x, ~, exitflag] = linprog(obj, D, d, N, b, lb, ub,linprog_options);
    case 2 % CPLEX      % exitflag ==  1: optimal, feasible  % exitflag == -2: infeasible
        [x, ~, exitflag] = cplexlp(obj, D, d, N, b, lb, ub);
    otherwise
        error('define a valid solver.');
end

if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
        || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
    feasible = 1;
    r = x;
else
    feasible = 0;
    return
end

if parsimonious == 1
fixed_fluxes = find(sum(logical(D),1));

% split network into forward and backward reactios
N2  = [N,-N];
ub2 = [ub;max([zeros(size(lb,1),1),-lb]')'];
lb2 = [max([zeros(size(lb,1),1),lb]')';-min([zeros(size(lb,1),1),ub]')'];
% fixate the already optimized fluxes
ub2(fixed_fluxes) = max([zeros(length(fixed_fluxes),1) r(fixed_fluxes)]')';
lb2(fixed_fluxes) = ub2(fixed_fluxes);
ub2(fixed_fluxes+cnap.numr) = -min([zeros(length(fixed_fluxes),1) r(fixed_fluxes)]')';
lb2(fixed_fluxes+cnap.numr) = ub2(fixed_fluxes+cnap.numr);
obj2 = ones(2*cnap.numr, 1);

switch solver
    case 0 % GLPK
        ctype = repmat('S',1,cnap.numis);
        [x2, ~, exitflag] =  glpk(  obj2,         N2, b, lb2, ub2, ctype);
    case 1 % linprog
        [x2, ~, exitflag] = linprog(obj2, [], [], N2, b, lb2, ub2, linprog_options);
    case 2 % CPLEX
        [x2, ~, exitflag] = cplexlp(obj2, [], [], N2, b, lb2, ub2);
    otherwise
        error('define a valid solver.');
end
if (exitflag == 1 && (solver == 2 || solver == 1)) ... CPLEX & Linprog
        || ((exitflag == 2 || exitflag == 5) && solver == 0) % GLPK
    r = x2(1:cnap.numr) - x2(cnap.numr+1:end);
else
    error('numerical error');
end
end
end