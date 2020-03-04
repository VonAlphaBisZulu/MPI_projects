function feasible = checkVectorFeasibility(N,lb,ub,V,v)
if nargin <= 3
    V = [];
    v = [];
end
    cgp= Cplex();
    cgp.Param.emphasis.numerical.Cur= 1;
    cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
    cgp.Model.A= N;
    cgp.Model.ub= ub;
    cgp.Model.lb= lb;
    cgp.Model.lhs= zeros(size(N,1), 1);
    cgp.Model.rhs= zeros(size(N,1), 1);
    cgp.Model.obj=zeros(size(N,2), 1);
    cgp.Model.sense= 'maximize';
    cgp.DisplayFunc=[];
    if ~isempty(v)
        cgp.addRows(-inf(size(V,1),1),V,v);
    end
    x= cgp.solve();
    if x.status == 3
        feasible = 0;
    else
        feasible = 1;
    end
end

