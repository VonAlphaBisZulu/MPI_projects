function valid = verify_mcs(cnap,modules,mcs,solver,verbose)
    if nargin < 5
        verbose = 0;
    end
    if nargin < 4
        solver = -1;
    end
    valid = nan(size(mcs,2),numel(modules));
    for i = 1:size(mcs,2)
        kos = find(mcs(:,i)==-1 | isnan(mcs(:,i)));
        A_kos = sparse( 1:2*numel(kos),[kos' kos'],...
            [ones(1,numel(kos)) -ones(1,numel(kos))],...
            2*numel(kos),cnap.numr);
        b_kos = zeros(2*numel(kos),1);
        for j = 1:numel(modules)
            A_ineq = [A_kos; modules{j}.V];
            b_ineq = [b_kos; modules{j}.v];
            switch modules{j}.type
                case 'lin_constraints'
                    fv = CNAoptimizeFlux(cnap,[],[],solver,verbose,0,A_ineq,b_ineq);
                    valid(i,j) = ~any(isnan(fv));
                case 'bilev_w_constr'
                    cnap.objFunc = modules{j}.c(:);
                    fv = CNAoptimizeFlux(cnap,[],[],solver,verbose,0,A_ineq,b_ineq);
                    A_ineq = [A_ineq; -modules{j}.c(:)'];
                    b_ineq = [b_ineq; -modules{j}.c(:)'*fv];
                    fv = CNAoptimizeFlux(cnap,[],[],solver,verbose,0,A_ineq,b_ineq);
                    valid(i,j) = ~any(isnan(fv));
                case 'yield_w_constr'
                case'raw'
            end
            switch modules{j}.sense
                case 'desired' % Desired Region
                case 'target' % Target Region
                    valid(i,j) = ~valid(i,j);
            end
        end
    end
end