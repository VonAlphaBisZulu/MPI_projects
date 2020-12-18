clear('cpx');
cpx.f     = cnap.objFunc;
cpx.sense = 'minimize';
cpx.Aineq = modules{1}.V;
cpx.bineq = modules{1}.v;
cpx.Aeq   = cnap.stoichMat;
cpx.beq   = zeros(cnap.nums,1);
cpx.lb    = cnap.reacMin;
cpx.ub    = cnap.reacMax;
cpx = Cplex(cpx);
cpx.Param.simplex.tolerances.feasibility.Cur = 0.00127291;
x = cpx.solve();
if x.status == 1
    disp(x.x(find(cnap.objFunc)));
    cnap.local.rb(:,1)=1:cnap.numr;
    cnap.local.rb(:,2)=x.x;
else
    disp(x.statusstring);
end
