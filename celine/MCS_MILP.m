classdef MCS_MILP
    % This class holds the MCS-MILP problem. The class MCS_solve is derived from
    % MCS_MILP and holds functions for the solution of the the MILP.
    % Classes like MCS_cplex or MCS_intlinprog are again derived from
    % MCS_solve and provide the solver interface.
    %
    % class constructor:
    % 
    %  obj = MCS_MILP(cnap,modules,koCost,kiCost,M,verbose)
    %
    % cnap        : a CellNetAnalyzer mass flow project
    % modules: cell-array of structs that define the desired and undesired flux 
    %          states in the form of (in)equality systems. Those systems describe
    %          subspaces of the original model (in "cnap") and can be forced to be 
    %          feasible or infeasible. There are three kinds of modules allowed:
    %           1. Inequalities ('lin_constraints'):
    %              e.g.: T r <= t.
    %           2. Optimality constraint with additional inequalities ('bilev_w_constr'):
    %              e.g.: c*r = optimal, T r <= t
    %           3. Yield constraints ('yield_w_constr'):
    %              e.g.: e*r/f*r <= t (Definition of V is not required, only v)
    %     fields:
    %        sense: 'desired' or 'target'
    %         type: 'lin_constraints', 'bilev_w_constr', 'yield_w_constr'
    %          V,v: Matrix/vector to specify additional linear constraints: V r <= v
    %               (e.g. T r <= t with 'target' or D r <= d with 'desired')
    %     module specific fields:
    %        lin_constraints: <none>
    %         bilev_w_constr: c: Inner optimization vector
    %         yield_w_constr: e: numerator of yield function,
    %                         f: denominator of yield function
    % koCost      : <double>[cnap.numr x 1] that indicates reactions that 
    %                can be knocked out. Notknockable reactions carry NaN. Knockable
    %                reactions carry a value that represents their knockout cost. By
    %                default all reactions are knockable.
    % kiCost      : <double>[cnap.numr x 1] that indicates reactions that 
    %                can be added or "knocked in". Not addable reactions carry NaN. 
    %                Addable reactions carry a value that represents their addition cost.
    %                If a reaction is set knockable and addable at the same time. It
    %                will be assumed that the reaction is addable. By default no
    %                reaction is addible.
    % M           : <double or bool> Specifies whether the big-M method or
    %                indicator constraints are used to link the decision
    %                variables z with the continuous part of the MILP.
    %                If M=0/false, indicator constraints are used. So 
    %                far this feature is only available with cplex
    %                and gurobi, but not with intlinprog or glpk. If
    %                M=true, MCS_MILP will choose a value. Other values
    %                of M may be used to directly specify the big-M value.
    % verbose     : Controls if this function reports to the command line.
    %
    % -----------------------
    % The MCS-MILP has the form:
    %
    %  max(c'x)
    %  A_ineq * x <= b_ineq
    %  A_eq   * x == b_eq
    %     lb <= x <= ub
    %  x(idx_z) % integer-variables in x (herein often named z)
    %
    % (optional) indicators (for details see function "link_z_indicators"):
    %  x(idx_z) == 1 -> A_indic * x <= b_indic
    %
    
    %% Issues/Discussion
    %
    %% 1. Find a good way to determine individual bigM for each constraint.
    %     - Problem should be tackled in link_z, where you call another function "best_M" to determine the
    %       ideal M for a particular constraint.
    %     - For now an universial M of e.g. 1000 is used
    %     - Issue 1: If M is too small, constraints (and thus reactions) cannot be comletely knocked out,
    %       hence some solutions are overlooked
    %     - e.g.: M = 1e2;
    %       -z*M + y1 + 3*y2 <= -2
    %       The solver tries to knock out: z=1
    %       y1 + 3*y2 <= -2 + z*M (1e2)
    %       y1 + 3*y2 <= 98
    %       -> if y1 and/or y2 have to take values in the order of magnitude 1e2 because
    %          of other constraints (like b'y <= -1), this constraint might again be binding at
    %          the higher value of 98 and prohibit a solution that might be beyond that point.
    %     - Issue 2: If M is too large, numerical issues occur, due to relaxation and tolerances.
    %       Either tolerances must then be set to very small values which slows down computation,
    %       or there are false positive solutions
    %     - e.g. integer tolerance is 1e-5 and M is 1e9:
    %       -z*M + y1 + 3*y2 <= -2
    %       The solver relaxes z=0 to z=1e-5:
    %       y1 + 3*y2 <= -2 + z*M (1e-5 * 1e9 = 1e4)
    %       y1 + 3*y2 <= 9998
    %       This can make the constraint non-biding without the need of setting z=1.
    %     - Idea: Find good M-values by optimization.
    %
    %% 2. Equilibrate problem
    %     - Preconditioning function was implemented. Improve / verify if it's helpful
    %
    %% 3. (Idea) Non-targetable z could be removed entirely from the problem.
    % That would reduce the number of integer variables. On the other hand, those
    % binary variables are already fixed to zero and might be removed by the solver
    % automatically. Not sure if worth the time.
    %
    %% 4. (Idea) Use z-variables as decision variables in an "either/or" way
    % z_map_vars, z_map_constr_eq, z_map_constr_ineq could be set up to make z-variables
    % work as knock-in indicators for some reactions and knock-out indicators
    % for others.
    %
    %% 5. (Idea) Maybe minimal sets of continuous interventions are possible
    % An intermediate variable v could be used as upper or lower bound for reactions
    % and z could then just indicate the number of non-zero v. Not sure how this
    % could work in target systems, because a smaller v in the primal translates
    % to a larger slack in the dual system.
    %
    
    properties (Access = 'public')
        A_ineq
        b_ineq
        A_eq
        b_eq
        c
        lb
        ub
        M % big M
        idx_z % indices of z-variables
        rownames_ineq
        rownames_eq
        colnames
        verbose % print messages
    end
    properties (Access = 'protected')
        indicators
        cost
        z_map_constr_ineq
        z_map_constr_eq
        z_map_vars
        z_inverted % inverted z (knock-ins)
        z_non_targetable % notknockable and non-addable z
        num_z
        cont_MILP
        cnap
    end
    properties (Access = 'private')
        num_modules
    end
    
    methods (Access = 'public')
        function obj = MCS_MILP(cnap,modules,koCost,kiCost,M,verbose)
            % Constructor of the MCS_MILP class:
            % Constructs the base problem: Initializes the constraint for the
            % sum of z-variables to be smaller than max Cost. Modules
            % are added on construction, but can also be added later.
            %
            % [-I        0        0        0        ]       <=   [-mincost  ]
            % [ I        0        0        0        ]       <=   [ maxcost  ]
            % [ H1_ineq  A1_ineq  0        0        ]  [ z] <=   [ b1_ineq  ]
            % [ H1_eq    A1_eq    0        0        ]  [x1]  =   [ b1_eq    ]
            % [ H2_ineq  0        A2_ineq  0        ]  [x2] <=   [ b2_ineq  ]
            % [ H2_eq    0        A2_eq    0        ]  [x3]  =   [ b2_eq    ]
            % [ H3_ineq  0        0        A3_ineq  ]       <=   [ b3_ineq  ]
            % [ H3_eq    0        0        A3_eq    ]        =   [ b3_eq    ]
            % [ ...                                 ] ...
            %
            % Construction parameters:
            % cnap: CNA project
            %
            % optional parameters:
            % modules: Strain design specification module. For detailed description, see "addModle"
            % koCost: Vector with knockout cost of each reaction. (nan for not-knockable)
            % kiCost: Vector with addition cost of each reaction.  (nan for non-addable)
            % maxCost: Maximal intervention costs (inf for no limit)
            % M: M used as big M for knocking out constraints or unbound reactions.
            % verbose: display messages
            %
            
            % process inputs
            n         = cnap.numr;
            obj.num_z = n ;
            obj.cnap  = CNAgetMFNetwork(cnap);
            if nargin < 2
                modules = [];
            end
            if nargin < 3 || isempty(kiCost)
                kiCost = nan(1,n);
            end
            if nargin < 4 || isempty(koCost)
                koCost = ones(1,n);
            end
            if nargin < 5 || isempty(M) || (islogical(M) && M)
                obj.M = 1e3;
            else
                obj.M = M;
            end
            if nargin < 6
                obj.verbose = 1;
            else
                obj.verbose = verbose;
            end
            
            % translate costs to inverted or non-targetable z-vectors
            koCost = koCost(:)'; % reshape vectors if necessary
            kiCost = kiCost(:)';
            koCost(~isnan(kiCost)) = nan;
            cost = koCost;
            cost(~isnan(kiCost)) = kiCost(~isnan(kiCost));
            obj.z_inverted = ~isnan(kiCost);      % these will be used later
            obj.z_non_targetable = isnan(cost);   % instead of koCost/kiCost
            % Set the sum of weighted z-variables as the upper limit and the
            % minimization of z-costs as the objective function
            cost(isnan(cost)) = 0;
            obj.cost = cost;
            obj.A_ineq = zeros(2,n);
            obj.A_ineq(1,:) = -cost; % min cost
            obj.A_ineq(2,:) =  cost; % max cost
            obj.b_ineq = [0; inf];
            obj.z_map_constr_ineq = zeros(n,2); % association between z and variables and inequalities
            obj.rownames_ineq = {'sum_z_min'; 'sum_z_max'};
            obj.colnames = cellfun(@(x) ['z' num2str(x)], num2cell(1:n), 'UniformOutput', false);
            obj.c = cost'; % objective function (minimize intervention costs)
            % for z lower bounds are zero, upper bounds are one (except for non-targetable reactions)
            obj.lb = zeros(n,1);
            obj.ub = ones(n,1);
            obj.ub(obj.z_non_targetable) = 0;
            obj.idx_z = 1:n;
            % Equality matrix is empty and will be extended by addModule when necessary
            obj.A_eq = nan(0,n);
            obj.b_eq = [];
            obj.z_map_constr_eq = zeros(n,0); % association between z and variables and equalities
            obj.rownames_eq = {};
            obj.num_modules = 0;
            obj.indicators = struct('A',zeros(0,obj.num_z),'b',[],'sense',[],'inverse',[],'type',[],'z',[],'epln',[]);
            % Initialize association between z and variables and variables
            obj.z_map_vars = zeros(n);
            for i = 1:numel(modules)
                obj = obj.addModule(modules{i},0);
            end
                        
            % save continous part of MILP for easy MCS validation
            cont_vars = setdiff(1:size(obj.A_ineq,2),obj.idx_z);
            obj.cont_MILP.A_ineq = obj.A_ineq(:,cont_vars);
            obj.cont_MILP.b_ineq = obj.b_ineq;
            obj.cont_MILP.A_eq   = obj.A_eq(:,cont_vars);
            obj.cont_MILP.b_eq   = obj.b_eq;
            obj.cont_MILP.lb     = obj.lb(cont_vars);
            obj.cont_MILP.ub     = obj.ub(cont_vars);
            obj.cont_MILP.c      = obj.c(cont_vars);
            obj.cont_MILP.z_map_constr_ineq = sparse(obj.z_map_constr_ineq);
            obj.cont_MILP.z_map_constr_eq   = sparse(obj.z_map_constr_eq);
            obj.cont_MILP.z_map_vars        = sparse(obj.z_map_vars(:,cont_vars));

            % 4. Link LP module to z-variables
            if obj.M == 0
                % Precondition matrix
%                 displ('Preconditioning MILP: finding variable bounds',obj.verbose);
%                    [ helper_lb,     helper_ub] = ...
%                 obj.bound_vars(  obj.A_ineq(3:end,cont_vars), obj.b_ineq(3:end), ...
%                                 obj.A_eq(:,cont_vars),       obj.b_eq, ...
%                                 obj.lb(cont_vars),     obj.ub(cont_vars), ...
%                                 obj.z_map_constr_ineq(:,3:end), obj.z_map_constr_eq,obj.z_map_vars(:,cont_vars));
%                 displ('Preconditioning MILP: scaling problem rows and columns',obj.verbose);
%                 [~, ~, ~,obj.A_ineq(3:end,cont_vars),  obj.b_ineq(3:end), ...
%                          obj.A_eq(:,cont_vars),        obj.b_eq, ...
%                          obj.lb(cont_vars),     obj.ub(cont_vars)] = ...
%                     preconditioning(obj,[], ...
%                          obj.A_ineq(3:end,cont_vars), obj.b_ineq(3:end), ...
%                          obj.A_eq(:,cont_vars),       obj.b_eq, ...
%                          obj.lb(cont_vars),     obj.ub(cont_vars), ...
%                          helper_lb,             helper_ub);
%                min(min(setdiff(abs([obj.A_ineq,obj.b_ineq;obj.A_eq, obj.b_eq;obj.lb',0; obj.ub',0]),[0,inf])))
%                max(max(setdiff(abs([obj.A_ineq,obj.b_ineq;obj.A_eq, obj.b_eq;obj.lb',0; obj.ub',0]),[0,inf])))
%                 preserve first two rows that define min and max cost, and don't scale columns of z variables to conserve bounds of +/- 1.
                obj = obj.link_z_indicators();
            else
                displ('Bounding Variables',obj.verbose);
                [ obj.lb(cont_vars),     obj.ub(cont_vars)] = ...
                                obj.bound_vars(  obj.A_ineq(3:end,cont_vars), obj.b_ineq(3:end), ...
                                obj.A_eq(:,cont_vars),       obj.b_eq, ...
                                obj.lb(cont_vars),     obj.ub(cont_vars), ...
                                obj.z_map_constr_ineq(:,3:end), obj.z_map_constr_eq,obj.z_map_vars(:,cont_vars));
                obj = obj.link_z();
            end      
        end
    end
    
    methods (Access = 'private')
        function obj = addModule(obj,module,precon)
            % Generating LP (A) and z-linking-matrix (H) for each module
            %
            % module: cell-array of structs that define the specific systems used in the final MILP.
            %         those system can be forced to be feasible or infeasible. There are three kinds
            %         of modules allowed:
            %           1. The wildtype model, constrainted with additional inequalities:
            %              e.g.: T r <= t.
            %           2. The wildtype model at a specific optimum, constrained with additional inequalities
            %              e.g.: c*r = optimal, T r <= t
            %           3. A yield range:
            %              e.g.: e*r/f*r <= t (Definition of V is not required, only v)
            %     fields:
            %        sense: 'desired' or 'target'
            %         type: 'lin_constraints', 'bilev_w_constr', 'yield_w_constr', 'raw'
            %          V,v: Matrix/vector to specify additional linear constraints: V r <= v
            %               (e.g. T r <= t with 'target' or D r <= d with 'desired')
            %     module specific fields:
            %        lin_constraints: <none>
            %         bilev_w_constr: c: Inner optimization vector
            %         yield_w_constr: e: numerator of yield function,
            %                         f: denominator of yield function
            %                    raw: A_ineq, b_ineq, A_eq, b_eq, lb, ub,
            %                         z_map_constr_ineq, z_map_constr_eq, z_map_vars
            %
            % precon: <true,false> if true, global MILP is preconditionned
            %
            
            obj.num_modules = obj.num_modules+1;
            z_map_constr_ineq_i = [];
            z_map_constr_eq_i = [];
            z_map_vars_i = [];
            if ~isfield(module,'lb') || isempty(module.lb)
                module.lb = obj.cnap.reacMin;
            end
            if ~isfield(module,'ub') || isempty(module.ub)
                module.ub = obj.cnap.reacMax;
            end
            if ~isfield(module,'V') || isempty(module.V)
                module.V = nan(0,obj.cnap.numr);
            end
            if ~isfield(module,'v') || isempty(module.v)
                module.v = [];
            end
            
            % 1. Construct LP for module
            switch module.type
                %% Classical MCS
                case 'lin_constraints'
                    [A_ineq_p, b_ineq_p, A_eq_p, b_eq_p, ~, lb_p, ub_p, z_map_constr_ineq_p, z_map_constr_eq_p, z_map_vars_p] = build_primal(obj,module.V,module.v,[],module.lb,module.ub);
                    %% Bilevel constraints
                case 'bilev_w_constr'
                    % 1. build primal w/ desired constraint (build_primal) - also store variable c
                    [A_ineq_v, b_ineq_v, A_eq_v, b_eq_v, c_v, lb_v, ub_v, z_map_constr_ineq_v, z_map_constr_eq_v, z_map_vars_v] = build_primal(obj,module.V,module.v,module.c,module.lb,module.ub);
                    % 2. build primal w/o desired constraint (build_primal) - store c_inner
                    [A_ineq_inner, b_ineq_inner, A_eq_inner, b_eq_inner, c_inner, lb_inner, ub_inner, z_map_constr_ineq_inner, z_map_constr_eq_inner, z_map_vars_inner] = build_primal(obj, [], [], module.c,module.lb,module.ub);
                    % 3. build dual from primal w/o desired constraint (build_dual w/o the farkas-option) - store c_inner_dual
                    [A_ineq_dual, b_ineq_dual, A_eq_dual, b_eq_dual, c_inner_dual, lb_dual, ub_dual, z_map_constr_ineq_dual, z_map_constr_eq_dual, z_map_vars_dual] = dualize(obj,A_ineq_inner, b_ineq_inner, A_eq_inner, b_eq_inner, c_inner, lb_inner, ub_inner,...
                        z_map_constr_ineq_inner ,z_map_constr_eq_inner, z_map_vars_inner);
                    % 4. connect primal w/ target region and dual w/o target region (i.e. biomass) via c = c_inner.
                    A_ineq_p=[A_ineq_v, zeros(size(A_ineq_v,1), size(A_ineq_dual,2)); zeros(size(A_ineq_dual,1), size(A_ineq_v,2)) A_ineq_dual];
                    b_ineq_p=[b_ineq_v; b_ineq_dual];
                    A_eq_p=[A_eq_v zeros(size(A_eq_v,1), size(A_eq_dual,2)); zeros(size(A_eq_dual,1), size(A_eq_v,2)) A_eq_dual; c_v, c_inner_dual];
                    b_eq_p=[b_eq_v; b_eq_dual;0];
                    lb_p=[lb_v; lb_dual];
                    ub_p=[ub_v; ub_dual];
                    
                    % 5. Update z-associations
                    z_map_vars_p=[z_map_vars_v z_map_vars_dual];
                    z_map_constr_ineq_p=[z_map_constr_ineq_v ,z_map_constr_ineq_dual];
                    z_map_constr_eq_p=[z_map_constr_eq_v, z_map_constr_eq_dual , zeros(obj.num_z,1)];
                    %% Yield constraints
                case 'yield_w_constr'
                    % numerator: e, denominator: f, yield: y
                    % change signs if denominator is negative
                    cnap1 = obj.cnap;
                    cnap1.objFunc =  module.f';
                    fv1 =  CNAoptimizeFlux(cnap1,[],[],-1,0,0,module.V,module.v);
                    cnap1.objFunc = -module.f';
                    fv2 = CNAoptimizeFlux(cnap1,[],[],-1,0,0,module.V,module.v);
                    den_min = module.f*fv1;
                    den_max = module.f*fv2;
                    if den_min * den_max >= 0 % check if denominator takes only positive or negative values and is ~=0
                        sign_denom = sign(den_min + den_max);
                        if sign_denom == 0
                            error('The divisor in a defined region seems to be (directly or indirectly) bound to zero. Please revise model and MCS setup.');
                        end
                    else
                        error('The divisor must not range positive AND negative values.');
                    end
                    module.e = module.e*sign_denom;
                    module.f = module.f*sign_denom;
                    
                    % Build LP
                    [A_ineq_v, b_ineq_v, A_eq_v, b_eq_v, ~, lb_v, ub_v, z_map_constr_ineq_v, z_map_constr_eq_v, z_map_vars_v] = build_primal(obj,module.V,module.v,[],module.lb,module.ub);
                    lb_nonzero = find(lb_v & ~isinf(lb_v));
                    lb_pos = lb_v >= 0;
                    ub_nonzero = find(ub_v & ~isinf(ub_v));
                    ub_neg = ub_v <= 0;
                    % prepare scalable upper and lower bounds
                    I_lb = full(sparse(1:length(lb_nonzero),lb_nonzero,-1,length(lb_nonzero),length(lb_v)));
                    I_ub = full(sparse(1:length(ub_nonzero),ub_nonzero, 1,length(ub_nonzero),length(ub_v)));
                    A_ineq_v = [[A_ineq_v; I_lb; I_ub], [-b_ineq_v; lb_v(lb_nonzero); -ub_v(ub_nonzero)]];
                    b_ineq_v = zeros(numel(b_ineq_v)+length(lb_nonzero)+length(ub_nonzero),1);
                    lb_v(lb_nonzero) = -inf;
                    lb_v(lb_pos) = 0;
                    ub_v(ub_nonzero) = inf;
                    ub_v(ub_neg) = 0;
                    lb_v(end+1) = 0;
                    ub_v(end+1) = inf;
                    lb_p = lb_v;
                    ub_p = ub_v;
                    % add constraint: denominator f equals 1 (or -1 if denominator is negative)
                    A_eq_p = [[A_eq_v; module.f] [-b_eq_v; 0]];
                    b_eq_p = [zeros(size(A_eq_v,1),1); sign_denom];
                    % add yield constraint: e*r/f*r <= y
                    % <=> (e-f*y)*r) <= 0
                    A_ineq_p = [A_ineq_v; (module.e*sign_denom - module.f*sign_denom*module.y), 0];
                    b_ineq_p = [b_ineq_v; 0];
                    
                    z_map_constr_ineq_p = [z_map_constr_ineq_v zeros(obj.num_z,length(lb_nonzero)+length(ub_nonzero)+1)];
                    z_map_constr_eq_p = [z_map_constr_eq_v zeros(obj.num_z,1)];
                    z_map_vars_p = [z_map_vars_v zeros(obj.num_z,1)];
                case'raw'
                    A_ineq_i = module.A_ineq;
                    b_ineq_i = module.b_ineq;
                    A_eq_i   = module.A_eq;
                    b_eq_i   = module.b_eq;
                    lb_i     = module.lb;
                    ub_i     = module.ub;
                    z_map_constr_ineq_i = module.z_map_constr_ineq;
                    z_map_constr_eq_i   = module.z_map_constr_eq;
                    z_map_vars_i        = module.z_map_vars;
                    if size(z_map_constr_ineq_i,1) ~= obj.num_z || size(z_map_constr_eq_i,1) ~= obj.num_z || size(z_map_vars_i,1) ~= obj.num_z
                        error('The number of z-variable must be identical to the number of reactions in the model');
                    end
                otherwise
                    error(['module type ' module.type ' unknown']);
            end
            % 2. Prepare module as target or desired
            switch module.sense
                case 'desired' % Desired Region
                    [A_ineq_i, b_ineq_i, A_eq_i, b_eq_i, lb_i, ub_i, z_map_constr_ineq_i, z_map_constr_eq_i] = reassign_lb_ub_from_ineq(obj,A_ineq_p, b_ineq_p, A_eq_p, b_eq_p, lb_p, ub_p, z_map_constr_ineq_p, z_map_constr_eq_p, z_map_vars_p);
                    z_map_vars_i        = z_map_vars_p;
                case 'target' % Target Region
                    c_p = zeros(1,size(A_eq_p,2));
                    [A_ineq_d, b_ineq_d, A_eq_d, b_eq_d, c_d, lb_i, ub_i, z_map_constr_ineq_d, z_map_constr_eq_d, z_map_vars_i] = dualize(obj,A_ineq_p, b_ineq_p, A_eq_p, b_eq_p, c_p, lb_p, ub_p, z_map_constr_ineq_p, z_map_constr_eq_p, z_map_vars_p);
                    [A_ineq_i, b_ineq_i, A_eq_i, b_eq_i, z_map_constr_ineq_i, z_map_constr_eq_i] = dual_2_farkas(obj,A_ineq_d, b_ineq_d,A_eq_d, b_eq_d, c_d, z_map_constr_ineq_d,z_map_constr_eq_d);
            end
            % 3. Add module to global MILP
            obj.z_map_constr_ineq   = [obj.z_map_constr_ineq,   z_map_constr_ineq_i];
            obj.z_map_constr_eq     = [obj.z_map_constr_eq,     z_map_constr_eq_i];
            obj.z_map_vars          = [obj.z_map_vars,          z_map_vars_i];
            obj.A_ineq  = [obj.A_ineq,                                  zeros(size(obj.A_ineq,1), size(A_ineq_i,2)) ; ...
                zeros(size(A_ineq_i,1), size(obj.A_ineq,2)), A_ineq_i];
            obj.b_ineq  = [obj.b_ineq;  b_ineq_i];
            obj.A_eq    = [obj.A_eq,                                zeros(size(obj.A_eq,1),size(A_eq_i,2)) ; ...
                zeros(size(A_eq_i,1), size(obj.A_eq,2)), A_eq_i];
            obj.b_eq    = [obj.b_eq;  b_eq_i];
            obj.c  = [obj.c;  zeros(size(A_ineq_i,2),1)];
            obj.lb = [obj.lb; lb_i];
            obj.ub = [obj.ub; ub_i];
            obj.rownames_ineq = [ obj.rownames_ineq ; cellfun(@(x) [num2str(obj.num_modules) '_mod_ineq_' num2str(x)], num2cell(1:size(A_ineq_i,1)), 'UniformOutput', false)'];
            obj.rownames_eq = [ obj.rownames_eq ; cellfun(@(x) [num2str(obj.num_modules) '_mod_eq_' num2str(x)], num2cell(1:size(A_eq_i,1)), 'UniformOutput', false)'];
            obj.colnames = [ obj.colnames , cellfun(@(x) [num2str(obj.num_modules) '_mod_var_' num2str(x)], num2cell(1:size(A_eq_i,2)), 'UniformOutput', false)];
        end
        
        function [A_ineq, b_ineq, A_eq, b_eq, c, lb, ub, z_map_constr_ineq, z_map_constr_eq, z_map_vars] = build_primal(obj,V,v,c,lb,ub)
            % Translates a CNA project and linear constraints to a problem of the
            %   standard form: A x = 0, A x <= b, lb <= x <= ub, min{c'x}.
            %   To only translate cnap to the standard form, set z_inverted to false and
            %   set z_non_targetable to true for all reactions;
            %
            % Output system: [A_ineq] . [x]' <= [b_ineq]
            %                [A_eq  ]         = [b_eq]
            %
            %                       lb_x <= x <= ub_x
            %
            % Input:
            % cnap: CNA project
            % V and v: Constraints of the form V r <= v
            %
            % Output:
            % z_map_vars: Matrix that connects the linear inequality systems to the global
            %    z-Variables, rows: z, associated with reactions/variables in the columns
            %    entries of "1" indicate that z is a knockout, entries of "-1" indicate
            %    that z=1 will be a "knock-in".
            % A_ineq, b_ineq, A_eq, b_eq, c: Linear Program of the form
            % A_ineq x <= b_ineq,   A_eq x = b_eq,   lb_x <= x < ub_x, min{x}.
            %
            if nargin < 2
                V = [];
                v = [];
            end
            if nargin < 4
                c  = obj.cnap.objFunc(:)';
            elseif isempty(c)
                c = zeros(1,obj.cnap.numr);
            end
            % fill matrices
            A_eq = obj.cnap.stoichMat;
            b_eq = zeros(obj.cnap.nums,1);
            if ~isempty(V) && ~isempty(v)
                A_ineq = V;
                b_ineq = v;
            else
                A_ineq = nan(0,obj.cnap.numr);
                b_ineq = nan(0,1);
            end
            if nargin < 5 || isempty(lb) || isempty(ub)
                lb = obj.cnap.reacMin;
                ub = obj.cnap.reacMax;
            end
            
            % map z-variables with continuous variables
            z_map_vars = eye(obj.num_z);
            z_map_vars(obj.z_non_targetable,:) = 0;
            z_map_vars(obj.z_inverted,:) = -z_map_vars(obj.z_inverted,:);
            z_map_constr_ineq = zeros(obj.num_z,size(A_ineq,1));
            z_map_constr_eq   = zeros(obj.num_z,size(A_eq,1));
            
            % Express positive lb or negative ub as inequalities to prevent
            % "non-convex" knockouts
            [A_ineq, b_ineq, lb, ub, z_map_constr_ineq] = prevent_boundary_knockouts(obj, A_ineq, b_ineq, lb, ub, z_map_constr_ineq, z_map_vars);
        end
        
        function [A_ineq, b_ineq, A_eq, b_eq, c, lb, ub, z_map_constr_ineq, z_map_constr_eq, z_map_vars] = dualize(obj,A_ineq_p, b_ineq_p, A_eq_p, b_eq_p, c_p, lb_p, ub_p,...
                z_map_constr_ineq_p, z_map_constr_eq_p, z_map_vars_p)
            % Translates a primal system to a dual system. The primal system must
            % be given in the standard form: A_ineq x <= b_ineq, A_eq x = b_eq, lb <= x < ub, min{x}.
            %
            % Variables translate to constraints:
            % x={R} ->   =
            % x>=0  ->  >= (new constraint is multiplied with -1 to translate to <=
            %               e.g. -A_i' y <= -c_i)
            % x<=0  ->  <=
            % Constraints translate to variables:
            % =     ->   y={R}
            % <=    ->   y>=0
            %
            %
            
            %
            % Consider that the following is not implemented:
            % In the case of (1) A x = b, (2) x={R}, (3) b~=0, Farkas' lemma is special,
            % because b'y ~= 0 is required to make the primal infeasible instead of b'y < 0.
            % 1. This does not occur very often.
            % 2. Splitting the equality into two inequalities that translate to y>=0
            %    would be posible, and yield b'y < 0 in the farkas' lemma.
            % Maybe splitting is required, but I actually don't think so. Using the
            % special case of b'y < 0 for b'y ~= 0 should be enough.
            
            if ~exist('z_map_vars_p','var') || isempty(z_map_vars_p)
                z_map_vars_p = zeros(obj.num_z,size(A_ineq_p,2));
            end
            if ~exist('z_map_constr_eq_p','var') || isempty(z_map_constr_eq_p)
                z_map_constr_eq_p = zeros(obj.num_z,size(A_eq_p,1));
            end
            if ~exist('z_map_constr_ineq_p','var') || isempty(z_map_constr_ineq_p)
                z_map_constr_ineq_p = zeros(obj.num_z,size(A_ineq_p,1));
            end
            
            % knockouts of variables and constraints must not overlap in the problem matrix
            if any(any(A_ineq_p(any(z_map_constr_ineq_p,1),any(z_map_vars_p,1)))) || any(any(A_eq_p(any(z_map_constr_eq_p,1),any(z_map_vars_p,1))))
                error('Knockable variables and knockable constraints overlap in the problem matrix. An error must have been made during the construction of the primal problem.');
            end
            
            n = obj.cnap.numr;
            if isempty(c_p)
                c_p = zeros(1,n);
            end
            
            % Translate inhomogenous bounds into inequality constraints
            lb_inh_bounds = lb_p~=0 & ~isinf(lb_p);
            ub_inh_bounds = ub_p~=0 & ~isinf(ub_p);
            x_geq0 = lb_p>=0 & ub_p >0;
            x_eR =   lb_p< 0 & ub_p >0;
            x_leq0 = lb_p< 0 & ub_p<=0;
            
            LB = full(sparse(1:sum(lb_inh_bounds),find(lb_inh_bounds),-1,sum(lb_inh_bounds),size(A_ineq_p,2)));
            UB = full(sparse(1:sum(ub_inh_bounds),find(ub_inh_bounds), 1,sum(ub_inh_bounds),size(A_ineq_p,2)));
            A_ineq_p = [ A_ineq_p; ...
                LB; ...
                UB ];
            b_ineq_p = [ b_ineq_p; ...
                -lb_p(lb_inh_bounds); ...
                ub_p(ub_inh_bounds)];
            
            % Translate into dual system
            A_ineq = [-A_eq_p(:,x_geq0)', -A_ineq_p(:,x_geq0)'; A_eq_p(:,x_leq0)', A_ineq_p(:,x_leq0)'];
            b_ineq = [c_p(x_geq0)' ; -c_p(x_leq0)'];
            A_eq = [A_eq_p(:,x_eR)', A_ineq_p(:,x_eR)'];
            b_eq = c_p(x_eR)';
            lb = [-inf(size(A_eq_p,1),1); zeros(size(A_ineq_p,1),1)];
            ub =   inf(size(A_eq_p,1)+size(A_ineq_p,1),1);
            c  = [b_eq_p; b_ineq_p]';
            
            % translate mapping of z-variables to rows instead of columns
            z_map_constr_ineq = [z_map_vars_p(:, x_geq0), z_map_vars_p(:, x_leq0)];
            z_map_constr_eq   = z_map_vars_p(:,x_eR);
            z_map_vars        = [z_map_constr_eq_p, z_map_constr_ineq_p, zeros(obj.num_z, size([LB;UB],1))];
            
            [A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_map_constr_ineq, z_map_constr_eq] = reassign_lb_ub_from_ineq(obj,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_map_constr_ineq, z_map_constr_eq, z_map_vars);
        end
        
        function [A_ineq, b_ineq, A_eq, b_eq, z_map_constr_ineq, z_map_constr_eq] = dual_2_farkas(~,A_ineq, b_ineq,A_eq, b_eq, c_dual, z_map_constr_ineq, z_map_constr_eq)
%             % add constraint b_prim'y or (c_dual'*y) <= -1;
%             A_ineq = [ A_ineq; c_dual];
%             b_ineq = [ b_ineq; -1];
%             z_map_constr_ineq = [z_map_constr_ineq, zeros(size(z_map_constr_ineq,1),1)];
            A_eq = [ A_eq; c_dual];
            b_eq = [ b_eq; -1];
            z_map_constr_eq = [z_map_constr_eq, zeros(size(z_map_constr_eq,1),1)];
            % it would also be possible (but ofc not necessary) to force (c_dual*y) == -1;
        end
        
        function obj = link_z_indicators(obj)
            % Linking z-variables to rows and columns of the module LP.
            % The LP (in A_ineq, b_ineq, A_eq, b_eq, lb, ub) is modified to
            % pair with the matrix G in a way that:
            %
            %  inequalities: z = 0 -> A_ineq_i x <= b_ineq_i
            %    equalities: z = 0 -> A_eq_i x   <= b_eq_i
            %     variables: z = 1 -> x_i = 0
            %
            % continu: lb_x <= x <= ub_x
            % integer: z = {0,1}
            %
            % The following matrices define the mapping between z variables and
            % constraints and variables (rows and columns) in the module LP
            % that are supposed to be knocked out or activated.
            %
            % z_map_constr_eq:   matrix z x rows(A_eq)
            % z_map_constr_ineq: matrix z x rows(A_ineq)
            % z_map_vars:        matrix z x columns(A_ineq or A_eq)
            %
            % Entries of  1 (e.g. in z_map_vars) mean that z=1 -> variable is knocked out
            % Entries of -1 (e.g. in z_map_vars) mean that z=0 -> variable is knocked out
            % Entries of  0 mark no relationship between a particular variable or constraint and the z-Variable
            %
            use_slack_vars = 0; % currently, indicator constraints are used without slack (v)
            knockable_ineq = any(obj.z_map_constr_ineq,1);
            knockable_eq   = any(obj.z_map_constr_eq,1);
            if use_slack_vars
                % add slack vars for knockable constraints and then link
                % slack variables to integer variables (z)
                % inequalities:
                obj.A_ineq(knockable_ineq,end+(1:sum(knockable_ineq)))  = -eye(sum(knockable_ineq));
                obj.A_eq(:,end+(1:sum(knockable_ineq))) = 0;
                obj.lb(end+(1:sum(knockable_ineq))) = 0;
                obj.ub(end+(1:sum(knockable_ineq))) = 1e4; % dummy bounds (could be replaced by bigM)
                obj.c(end+(1:sum(knockable_ineq))) = 0;
                obj.z_map_vars = [obj.z_map_vars -obj.z_map_constr_ineq(:,knockable_ineq)];
                obj.colnames(end+(1:sum(knockable_ineq))) = cellstr(strcat('slack_ineq_',strtrim(cellstr(num2str((1:sum(knockable_ineq))')))))';
                % Indicator constraints
                % equalities:
                obj.A_ineq(:,end+(1:sum(knockable_eq))) = 0;
                obj.A_eq(knockable_eq,end+(1:sum(knockable_eq))) = -eye(sum(knockable_eq));
                obj.lb(end+(1:sum(knockable_eq))) = -1e4;
                obj.ub(end+(1:sum(knockable_eq))) = 1e4;
                obj.c(end+(1:sum(knockable_eq))) = 0;
                obj.z_map_vars = [obj.z_map_vars -obj.z_map_constr_eq(:,knockable_eq)];
                obj.colnames(end+(1:sum(knockable_eq))) = cellstr(strcat('slack_eq_',strtrim(cellstr(num2str((1:sum(knockable_eq))')))))';
            else
                % Use indicator constraints to directly control the
                % activity of constraints via z.
                % inequalities:
                indicators_ineq.name    = strcat(obj.rownames_ineq(knockable_ineq),'_z');
                indicators_ineq.A       = obj.A_ineq(knockable_ineq,:);
                indicators_ineq.b       = obj.b_ineq(knockable_ineq);
                indicators_ineq.sense   = repmat('L',sum(knockable_ineq),1);
                indicators_ineq.inverse = double(sum(obj.z_map_constr_ineq(:,knockable_ineq),1) > 0)';
                indicators_ineq.type    = 1*ones(sum(knockable_ineq),1);
                indicators_ineq.epln    = 0.0*ones(sum(knockable_ineq),1);
                [z,~] = find(obj.z_map_constr_ineq(:,knockable_ineq));
                indicators_ineq.z       = z;
                % remove z-associated constraints from matrix, because they're now be treated by indicator constraints
                obj.A_ineq = obj.A_ineq(~knockable_ineq,:);
                obj.b_ineq = obj.b_ineq(~knockable_ineq);
                obj.rownames_ineq = obj.rownames_ineq(~knockable_ineq);

                % equalities
                indicators_eq.name    = strcat(obj.rownames_eq(knockable_eq),'_z');
                indicators_eq.A       = obj.A_eq(knockable_eq,:);
                indicators_eq.b       = obj.b_eq(knockable_eq);
                indicators_eq.sense   = repmat('E',sum(knockable_eq),1);
                indicators_eq.inverse = double(sum(obj.z_map_constr_eq(:,knockable_eq),1) > 0)';
                indicators_eq.type    = 1*ones(sum(knockable_eq),1);
                indicators_eq.epln    = 0.0*ones(sum(knockable_eq),1);
                [z,~] = find(obj.z_map_constr_eq(:,knockable_eq));
                indicators_eq.z       = z;
                % remove constraints from matrix, because they're now treated by indicator constraints
                obj.A_eq = obj.A_eq(~knockable_eq,:);
                obj.b_eq = obj.b_eq(~knockable_eq);
                obj.rownames_eq = obj.rownames_eq(~knockable_eq);
            end
            
            inverted      = any(obj.z_map_vars<0,1);
            var_lb_inf    = any(obj.z_map_vars,1) & obj.lb' ~= 0 & isinf(obj.lb)' & obj.ub' == 0;
            var_ub_inf    = any(obj.z_map_vars,1) & obj.ub' ~= 0 & isinf(obj.ub)' & obj.lb' == 0;
            var_lb_ub_inf = any(obj.z_map_vars,1)  & isinf(obj.ub)' & isinf(obj.lb)';
            
            % unbound variables
            knockable_vars          = any(var_lb_inf | var_ub_inf | var_lb_ub_inf,1);
            indicators_vars.name    = strcat(obj.colnames(knockable_vars)','_z');
            indicators_vars.A       = full(sparse(1:sum(knockable_vars),find(knockable_vars),1,sum(knockable_vars),numel(knockable_vars)));
            indicators_vars.b       = zeros(sum(knockable_vars),1);
            indicators_vars.sense   = repmat(' ',sum(knockable_vars),1);
            indicators_vars.sense(var_lb_ub_inf(knockable_vars)) = 'E';
            indicators_vars.sense(var_ub_inf(knockable_vars))    = 'L';
            indicators_vars.sense(var_lb_inf(knockable_vars))    = 'G';
            indicators_vars.inverse = double(inverted(var_lb_inf | var_ub_inf | var_lb_ub_inf))';
            indicators_vars.type    = 1*ones(sum(knockable_vars),1);
            indicators_vars.epln    = 0.0*ones(sum(knockable_vars),1);
            [z,~] = find(obj.z_map_vars(:,var_lb_inf | var_ub_inf | var_lb_ub_inf));
            indicators_vars.z       = z;
            
            % all indicator constraints pooled
            if use_slack_vars
                obj.indicators.A      =  indicators_vars.A      ;
                obj.indicators.b      =  indicators_vars.b      ;
                obj.indicators.sense  =  indicators_vars.sense  ;
                obj.indicators.inverse=  indicators_vars.inverse;
                obj.indicators.type   =  indicators_vars.type   ;
                obj.indicators.epln   =  indicators_vars.epln   ;
                obj.indicators.z      =  indicators_vars.z      ;
                obj.indicators.name   =  indicators_vars.name   ;
            else
                obj.indicators.A      =  [indicators_ineq.A;         indicators_eq.A;        indicators_vars.A      ];
                obj.indicators.b      =  [indicators_ineq.b;         indicators_eq.b;        indicators_vars.b      ];
                obj.indicators.sense  =  [indicators_ineq.sense;     indicators_eq.sense;    indicators_vars.sense  ];
                obj.indicators.inverse=  [indicators_ineq.inverse;   indicators_eq.inverse;  indicators_vars.inverse];
                obj.indicators.type   =  [indicators_ineq.type;      indicators_eq.type;     indicators_vars.type   ];
                obj.indicators.epln   =  [indicators_ineq.epln;      indicators_eq.epln;     indicators_vars.epln   ];
                obj.indicators.z      =  [indicators_ineq.z;         indicators_eq.z;        indicators_vars.z      ];
                obj.indicators.name   =  [indicators_ineq.name;      indicators_eq.name;     indicators_vars.name   ];
            end
            % bound variables are associated with z directly
            var_lb_M = any(obj.z_map_vars,1) & obj.lb' ~= 0 & ~isinf(obj.lb)'; % mapping is only required if the lower bound isn't already zero (inf is 
            var_ub_M = any(obj.z_map_vars,1) & obj.ub' ~= 0 & ~isinf(obj.ub)'; % already covered by indicator constraints) same for the upper bound.
            LB = full(sparse(1:sum(var_lb_M),find(var_lb_M),-1,sum(var_lb_M),size(obj.A_ineq,2))); % continuous part
            LB(:,1:obj.num_z) = -obj.z_map_vars(:,var_lb_M)'.*obj.lb(var_lb_M); % bounds are multiplied with integer variables
            UB = full(sparse(1:sum(var_ub_M),find(var_ub_M), 1,sum(var_ub_M),size(obj.A_ineq,2))); % continuous part
            UB(:,1:obj.num_z) = obj.z_map_vars(:,var_ub_M)'.*obj.ub(var_ub_M); % bounds are multiplied with integer variables
            obj.A_ineq = [obj.A_ineq; LB; UB];
            obj.b_ineq = [obj.b_ineq; -obj.lb(var_lb_M).*~inverted(var_lb_M)'; obj.ub(var_ub_M).*~inverted(var_ub_M)'];
            obj.rownames_ineq 	= [ obj.rownames_ineq ; strcat(obj.colnames(var_lb_M)','_z_lb'); strcat(obj.colnames(var_ub_M)','_z_ub')];
        end
        function obj = link_z(obj)
            % Linking z-variables to rows and columns of the module LP.
            % The LP (in A_ineq, b_ineq, A_eq, b_eq, lb, ub) is modified to
            % pair with the matrix G in a way that:
            %
            %  [G_ineq   A_ineq] . [z]' <= [b_ineq]
            %  [0        A_eq  ]   [x]   = [b_eq]
            %
            % continu: lb_x <= x <= ub_x
            % integer: z = {0,1}
            %
            % The following matrices define the mapping between z variables and
            % constraints and variables (rows and columns) in the module LP
            % that are supposed to be knocked out or activated.
            %
            % z_map_constr_eq:   matrix z x rows(A_eq)
            % z_map_constr_ineq: matrix z x rows(A_ineq)
            % z_map_vars:        matrix z x columns(A_ineq or A_eq)
            %
            % Entries of  1 (e.g. in z_map_vars) mean that z=1 -> variable is knocked out
            % Entries of -1 (e.g. in z_map_vars) mean that z=0 -> variable is knocked out
            % Entries of  0 mark no relationship between a particular variable or constraint and the z-Variable
            %
            
            z_map_constr_ineq_i = -obj.z_map_constr_ineq; % row knockouts must be treated inversely to variable knockouts
            z_map_constr_eq_i = -obj.z_map_constr_eq; % row knockouts must be treated inversely to variable knockouts
            % Implement constraint-knockouts: When equality constraints are mapped to z, they
            % must first be translated to two inequalities.
            eq2ineq = any(z_map_constr_eq_i,1);
            z_map_constr_ineq_i = [z_map_constr_ineq_i, z_map_constr_eq_i(:,eq2ineq), z_map_constr_eq_i(:,eq2ineq)];
            A_ineq_z = [ obj.A_ineq ; ...
                -obj.A_eq(eq2ineq,:); ... knockable equalities must be translated to
                obj.A_eq(eq2ineq,:)]; %  inequalities.
            b_ineq_z = [ obj.b_ineq; ...
                -obj.b_eq(eq2ineq) ; ...
                obj.b_eq(eq2ineq)];
            obj.A_eq = obj.A_eq(~eq2ineq,:);
            obj.b_eq = obj.b_eq(~eq2ineq); % After this, z_map_constr_eq is integrated into z_map_constr_ineq and no longer required
            
            obj.rownames_ineq 	= [ obj.rownames_ineq ; strcat(obj.rownames_eq(eq2ineq),'_ge'); strcat(obj.rownames_eq(eq2ineq),'_le')];
            obj.rownames_eq     = obj.rownames_eq(~eq2ineq);
            
            
            % Implement variable-knockouts: Translate them to constraint-knockouts
            var_ko_lb = any(obj.z_map_vars,1) & obj.lb' ~= 0; % mapping is only required if the lower bound isn't already zero
            var_ko_ub = any(obj.z_map_vars,1) & obj.ub' ~= 0; % same for the upper bound
            LB = full(sparse(1:sum(var_ko_lb),find(var_ko_lb),-1,sum(var_ko_lb),size(A_ineq_z,2)));
            UB = full(sparse(1:sum(var_ko_ub),find(var_ko_ub), 1,sum(var_ko_ub),size(A_ineq_z,2)));
            A_ineq_z = [A_ineq_z ; LB; UB]; % prepare the contiuous part of the knockouts (quasi-identity-matrix)
            b_ineq_z  = [b_ineq_z; zeros(sum(var_ko_lb)+sum(var_ko_ub),1)];
            % variable knockouts are now also mapped to constraint knockouts.
            z_map_constr_ineq_i = [z_map_constr_ineq_i, obj.z_map_vars(:,var_ko_lb), obj.z_map_vars(:,var_ko_ub)];
            obj.rownames_ineq 	= [ obj.rownames_ineq ; strcat(obj.colnames(var_ko_lb)','_z_lb'); strcat(obj.colnames(var_ko_ub)','_z_ub')];
            
            M_ineq = zeros(size(A_ineq_z,1),1);
            idx_cont = setdiff(1:size(obj.A_eq,2),obj.idx_z);
            
            if exist('computeM.m','file')
                computeM = evalin('base','@computeM'); % If external computeM function exists, use that one instead of the built-in
            else
                computeM = evalin('caller','@computeM'); % Compute bigM - constraint-wise if possible.
            end
            [M_ineq(3:end), b_eq0] = computeM(obj,A_ineq_z(3:end,idx_cont), b_ineq_z(3:end,:), obj.A_eq(:,idx_cont), obj.b_eq, obj.lb(idx_cont), obj.ub(idx_cont), z_map_constr_ineq_i(:,3:end));
            obj.b_eq = b_eq0;
            obj.A_ineq = A_ineq_z;
            
            if any(sum(abs(z_map_constr_ineq_i),1)>1)
                error('Every variable or consraint can only be linked to one z-vector.');
            end
            
            % use M as matrix coefficients multiplied with z and for rhs if necessary
            z_coeff = M_ineq; % coefficients in the z-Matrix (e.g. ub, lb, -M, M (that represent -inf, inf))
            b_ineq_add = M_ineq; % values that are added to b_ineq (can also be -M, M etc.)
            
            % identify inverted knockouts and change coefficients for z, and terms added to b
            inverted = any(z_map_constr_ineq_i<0,1);
            b_ineq_add(inverted) = 0;
            z_coeff(inverted) = -z_coeff(inverted);
            z_map_constr_ineq_i = abs(z_map_constr_ineq_i);
            
            % build the z-associating part of the matrix A_ineq (in the description of the constructor called H) to complete the row knockout linking
            H = z_map_constr_ineq_i'.*z_coeff;
            obj.A_ineq(any(z_map_constr_ineq_i,1),1:obj.num_z) = H(any(z_map_constr_ineq_i,1),:); % for knockable reactions, multiply integer variables with bounds
            obj.b_ineq = b_ineq_z+b_ineq_add;
        end
        
        function [lb_new, ub_new] = bound_vars(obj,A_ineq,b_ineq,A_eq,b_eq,lb,ub,z_map_constr_ineq,z_map_constr_eq,z_map_vars)
            % This function identifies (implicit) bounds on all variables.
            % Those bounds are accounted for in preconditioning to make
            % sure that variables don't take unnecessarily small values.
            lp.A_ineq = A_ineq(~any(z_map_constr_ineq,1),:);
            lp.b_ineq = b_ineq(~any(z_map_constr_ineq,1));
            lp.A_eq   = A_eq(~any(z_map_constr_eq,1),:);
            lp.b_eq   = b_eq(~any(z_map_constr_eq,1));
            lp.lb     = lb;
            lp.ub     = ub;
            lp.c      = zeros(size(A_ineq,2),1);
            if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
                   (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && ~isempty(gcp('nocreate')) %#ok<DCRENAME>
                numworkers = getfield(gcp('nocreate'),'NumWorkers');
            else
                numworkers = 0;
            end
            lb_new = lp.lb;
            ub_new = lp.ub;
            solveLP = @obj.solveLP;
            parfor (i = 1:numel(ub),numworkers)
                lp_1 = lp;
                if lb_new(i) ~= 0
                    lp_1.c(i) = 1;
                    x = solveLP(lp_1);
                    if ~isempty(x)
                        lb_new(i) = x(i);
                    end
                end
                if lb_new(i) ~= 0
                    lp_1.c(i) = -1;
                    x = solveLP(lp_1);
                    if ~isempty(x)
                        ub_new(i) = x(i);
                    end
                end
            end
        end
        
        function [c, factors_rows, factors_columns, A_ineq, b_ineq, A_eq, b_eq, lb, ub] = preconditioning(obj, c, A_ineq, b_ineq, A_eq, b_eq, lb, ub, lb_help, ub_help, fixed_ineq, fixed_eq, fixed_cols)
            % LP preconditioning to avoid large orders of magnitude between the matrix
            % entries. Matrix columns and rows are scaled with individual factors.
            % The scaling is optimizes all entries of the final problem to be as large
            % as possible, but smaller equal one. Furthermore scaling factors for rows and
            % columns that should be exempt from this optimization can be defined.
            %
            % Uses logarithm to transform qp into lp
            %
            % 1e-6: row/column factors are rounded to 6 decimals. inf: no rounding
            precision = inf; % define the number of decimals for row- and column factors
            min_oom = 0; % define here smallest order of magnitude
            if isempty(c)
                c = zeros(size(A_ineq,2),1);
            end
            if nargin < 7
                lb = [];
                ub = [];
            end
            if nargin < 9
                lb_help = lb;
                ub_help = ub;
            end
            if nargin < 11
                fixed_ineq = [];
            end
            if nargin < 12
                fixed_eq = [];
            end       
            if nargin < 13
                fixed_cols = [];
            end
            num_x = size(A_ineq,2);
            range_ineq  =   (1:size(b_ineq,1));
            range_eq    =      size(b_ineq,1) + (1:size(b_eq,1));
            range_lb    =      size(b_ineq,1) +    size(b_eq,1) + 1:size(A_ineq,2);
            range_ub    =      size(b_ineq,1) +    size(b_eq,1) +   size(A_ineq,2) + 1:size(A_ineq,2);
            % build matrix of entire LP
            A = [ A_ineq, b_ineq; A_eq, b_eq; -eye(num_x), -lb_help; eye(num_x), ub_help];
            fixed_rows = [fixed_ineq, size(A_ineq,1)+fixed_eq, range_lb, range_ub];
            % prepare preconditioning matrix
            num_rows = size(A,1);
            num_cols = size(A,2);
            % log10 the matrix
            A_log = log10(abs(A));
            % remove non-scalable values from the matrix
            A_log(fixed_rows,fixed_cols) = -inf;
            % build LP
            % only include entries that can be scaled and were not 0 or inf (-inf or inf after log):
            [rows, cols] = find(~isinf(A_log)); % find these entries
            num_entries  = length(rows);
            cols = cols+num_rows;
            fixed_cols = fixed_cols+num_rows;
            % Translate quadratic problem into LP via logarithm
            % g^-1 <= a_ij*fc*fr <= g
            % after logarithm:
            % -ln_g <= ln_a + ln_fc + ln_fr <=  ln_g
            % -> (1)  -ln_g - ln_fc - ln_fr <=  ln_a
            % -> (2)  -ln_g + ln_fc + ln_fr <= -ln_a
            % this is the right hand side:
            b_log = [A_log(~isinf(A_log)); -A_log(~isinf(A_log))];
            % this is a list of the "coordinates" of each scalable entry in A
            ref = [rows cols; rows cols]';
            % In the preconditioning LP, these values are added twice. Once
            % for the upper bound, once for the lower bound problem.
            entries = [-ones(2*num_entries,1); ones(2*num_entries,1)];
            A_log = sparse( repelem(1:2*num_entries,1,2),...
                ref(:),...
                entries);
            % setup lb, ub and objective
            c_log = zeros(size(A_log,2),1);
            lb_log = -inf(size(c_log));
            ub_log =  inf(size(c_log));
            lb_log([fixed_rows fixed_cols]) = 0;
            ub_log([fixed_rows fixed_cols]) = 0;
            % Solve to find largest global lower bound ln_g (and smallest upper bound ln_h)
            % The difference between g and h will determine the largest
            % order-of-magnitude-difference in the lp.
            A_log1  = [ A_log, ...
                      [-ones(num_entries,1); -ones(num_entries,1)]]; % ln_g
            lb_log1 = [lb_log; 0  ];
            ub_log1 = [ub_log; inf];
            c_log1  = [c_log;  1  ];
            % Options
            x = obj.solveLP(struct('A_ineq',A_log1,'b_ineq',b_log,'A_eq',[],'b_eq',[],'lb',lb_log1,'ub',ub_log1,'c',c_log1));
%             disp(x(end)) % show lower bound
            % Solve again to shift all matrix entries as close as possible to 1
            % by minimizing | ln_a + ln_fc + ln_fr |
            % trough: min_oom <= l <=  ln_a + ln_fc + ln_fr <= h % min_oom is the lower bound on the order of magnitude
            % -> (1)  l - ln_fc - ln_fr <=  ln_a
            % -> (2) -h + ln_fc + ln_fr <= -ln_a
            % -> max(min_oom,-k) <= l <= 0
            % -> 0 <= h <= max(2*k-oom_min, k)
            A_log2  = [ A_log, [ eye(num_entries) ; zeros(num_entries) ], [ zeros(num_entries) ; -eye(num_entries) ]];
            % bound k (via k <= ln_h) with the global bound ln_g found in the last step
            lb_log2 = [lb_log; max(min_oom,-x(end))*ones(num_entries,1); zeros(num_entries,1) ];
            ub_log2 = [ub_log; zeros(num_entries,1);                      max(2*x(end)+min_oom,x(end))*ones(num_entries,1)];
            % minimize the sum of all k (maybe minimizing the sum of squares would be better?)
            c_log2  = [c_log;-ones(num_entries,1);ones(num_entries,1)];
            x = obj.solveLP(struct('A_ineq',A_log2,'b_ineq',b_log,'A_eq',[],'b_eq',[],'lb',lb_log2,'ub',ub_log2,'c',c_log2));
            % translate factors back
            factors_rows    = 10.^(x(1:num_rows));
            factors_columns = 10.^(x(num_rows+(1:num_cols)));
            % round
            if ~isinf(precision)
                factors_rows    = round(factors_rows,precision,'significant'); % floor after 5 decimals
                factors_columns = round(factors_columns,precision,'significant');
            end
            A = (A'.*factors_columns)'.*factors_rows; % multiply rows and columns
            
            A_ineq = (A_ineq'.*factors_columns(1:end-1))'.*factors_rows(range_ineq);
            b_ineq = b_ineq.*factors_columns(end).*factors_rows(range_ineq);
            A_eq   = (A_eq'.*factors_columns(1:end-1))'.*factors_rows(range_eq);
            b_eq   = b_eq.*factors_columns(end).*factors_rows(range_eq);
            lb     = lb.*factors_columns(end)./factors_columns(1:end-1);
            ub     = ub.*factors_columns(end)./factors_columns(1:end-1);
            
            % map back to input
%             A_ineq = A(1:size(A_ineq,1),1:end-1);
%             b_ineq = A(1:size(b_ineq,1),end);
%             A_eq   = A(size(A_ineq,1)+(1:size(A_eq,1)),1:end-1);
%             b_eq   = A(size(A_ineq,1)+(1:size(A_eq,1)),end);
%             lb_A   = diag(A(size([A_ineq;A_eq],1)+(1:num_x),1:end-1));
%             lb_b   = A(size([A_ineq;A_eq],1)+(1:num_x),end);
%             lb     = lb_b./lb_A;
%             ub_A   = diag(A(size([A_ineq;A_eq],1)+num_x+(1:num_x),1:end-1));
%             ub_b   = A(size([A_ineq;A_eq],1)+num_x+(1:num_x),end);
%             ub     = ub_b./ub_A;
            % adapt objective function
            c = c./factors_columns(1:end-1);
        end
        
        function [A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_map_constr_ineq, z_map_constr_eq] = reassign_lb_ub_from_ineq(~,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_map_constr_ineq, z_map_constr_eq, z_map_vars)
            % this function searches for bounds on variables (single entries in A_ineq or A_eq). Those
            % constraints are removed and instead put into lb and ub. This
            % is useful to avoid unnecessary bigM constraints on variable
            % bounds when true upper and lower bounds exist instead.
            % To avoid interference with the knock-out
            % logic, negative ub and positive ub are not translated.
            n = size(A_ineq,2);
            lb = num2cell(lb);
            ub = num2cell(ub);
            % translate entries to lb or ub
            % find all entries in A_ineq and filter for rows with only one entry
            [row_ineq,  ~] = find(A_ineq);
            [num_occ_ineq, row_ineq] = hist(row_ineq,unique(row_ineq));
            %                                                                don't include positive lb or negative ub
            var_bound_constraint_ineq = setdiff(row_ineq(num_occ_ineq == 1)',find(b_ineq<0));
            % don't include knockable constraints
            var_bound_constraint_ineq = setdiff(var_bound_constraint_ineq, find(any(z_map_constr_ineq,1)));
            for i = var_bound_constraint_ineq(:)'
                idx_r = find(A_ineq(i,:));
                if A_ineq(i,idx_r)>0 % upper bound constraint
                    ub{idx_r} = [ub{idx_r} b_ineq(i)/A_ineq(i,idx_r)];
                else % lower bound constraint
                    lb{idx_r} = [lb{idx_r} b_ineq(i)/A_ineq(i,idx_r)];
                end
            end
            % find all entries in A_eq and filter for rows with only one entry
            [row_eq  ,  ~] = find(A_eq);
            [num_occ_eq,   row_eq]  = hist(row_eq,unique(row_eq));
            var_bound_constraint_eq = row_eq(num_occ_eq == 1)';
            % don't include knockable constraints
            var_bound_constraint_eq = setdiff(var_bound_constraint_eq, find(any(z_map_constr_eq,1)));
            % Also partly set lb or ub derived from equality constraints, for instance:
            % If x =  5, set ub = 5 and keep the inequality constraint -x <= -5.
            % If x = -5, set lb =-5 and keep the inequality constraint  x <=  5.
            A_ineq_new = [];
            b_ineq_new = [];
            for i = var_bound_constraint_eq(:)'
                idx_r = find(A_eq(i,:));
                if any(z_map_vars(:,idx_r)) % if knockable
                    if A_eq(i,idx_r) * b_eq(i) >0 % upper bound constraint
                        ub{idx_r} = [ub{idx_r} b_eq(i)/A_eq(i,idx_r)];
                        A_ineq_new = [A_ineq_new; -A_eq(i,:)];
                        b_ineq_new = [b_ineq_new; -b_eq(i)];
                    elseif A_eq(i,idx_r) * b_eq(i) <0 % lower bound constraint
                        lb{idx_r} = [lb{idx_r} b_eq(i)/A_eq(i,idx_r)];
                        A_ineq_new = [A_ineq_new;  A_eq(i,:)];
                        b_ineq_new = [b_ineq_new;  b_eq(i)];
                    else
                        lb{idx_r} = [lb{idx_r} 0];
                        ub{idx_r} = [ub{idx_r} 0];
                    end
                else % if notknockable
                    lb{idx_r} = [lb{idx_r} b_eq(i)/A_eq(i,idx_r)];
                    ub{idx_r} = [ub{idx_r} b_eq(i)/A_eq(i,idx_r)];
                end
            end
            for i = 1:n
                if any(z_map_vars(:,i))   % if knockable, select losest bounds
                    ub{i} = max(ub{i}(~isinf(ub{i}))); % select largest non-inf upper bound
                    lb{i} = min(lb{i}(~isinf(lb{i}))); % select smallest non-inf lower bound
                else % if not knockable, select tightest bounds
                    ub{i} = min(ub{i}(~isinf(ub{i}))); % select largest non-inf upper bound
                    lb{i} = max(lb{i}(~isinf(lb{i}))); % select smallest non-inf lower bound
                end
            end
            ub(cellfun(@isempty,ub)) = { inf};
            lb(cellfun(@isempty,lb)) = {-inf};
            ub = cell2mat(ub);
            lb = cell2mat(lb);
            
            if any(lb>ub)
                error('Variables have contradictory bounds/constraints.');
            end
            % remove constraints that became redundant
            A_ineq(var_bound_constraint_ineq,:)            = [];
            b_ineq(var_bound_constraint_ineq)              = [];
            z_map_constr_ineq(:,var_bound_constraint_ineq) = [];
            A_eq(var_bound_constraint_eq,:)                = [];
            b_eq(var_bound_constraint_eq)                  = [];
            z_map_constr_eq(:,var_bound_constraint_eq)     = [];
            % add equality constraints that transformed to inequality constraints
            A_ineq            = [A_ineq; A_ineq_new];
            b_ineq            = [b_ineq; b_ineq_new];
            z_map_constr_ineq = [z_map_constr_ineq, zeros(size(z_map_constr_ineq,1),size(A_ineq_new,1))];
        end
        function [A_ineq_i, b_ineq_i, lb_i, ub_i, z_map_constr_ineq_i] = prevent_boundary_knockouts(obj, A_ineq, b_ineq, lb, ub, z_map_constr_ineq, z_map_vars)
            lb_i = lb;
            ub_i = ub;
            A_ineq_i = A_ineq;
            b_ineq_i = b_ineq;
            z_map_constr_ineq_i = z_map_constr_ineq;
            for i = 1:size(A_ineq,2)
                if any(z_map_vars(:,i)==1,1) && lb(i)>0
                    A_ineq_i = [A_ineq_i; full(sparse(1,i,-1,1,size(A_ineq,2)))];
                    b_ineq_i = [b_ineq_i; -lb(i)];
                    z_map_constr_ineq_i = [z_map_constr_ineq_i, zeros(obj.num_z,1)];
                    lb_i(i) = 0;
                end
                if any(z_map_vars(:,i)==1,1) && ub(i)<0
                    A_ineq_i = [A_ineq_i; full(sparse(1,i, 1,1,size(A_ineq,2)))];
                    b_ineq_i = [b_ineq_i; ub(i)];
                    z_map_constr_ineq_i = [z_map_constr_ineq_i, zeros(obj.num_z,1)];
                    ub_i(i) = 0;
                end
            end
        end
        function [M_ineq, b_eq0] = computeM(obj,A_ineq, b_ineq, A_eq, b_eq, lb, ub, z_ineq)
            b_eq=b_eq0;
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
            displ('Trying to find big-M values through optimization (MCS_MILP built-in function)',obj.verbose);
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
                    M_ineq(i) = -lp.c*x; % maybe + epsilon?
                else
                    M_ineq(i) = inf;
                end
            end
            M_ineq(isinf(M_ineq)) = obj.M;
            M_ineq(isnan(M_ineq)) = 0;
        end
    end
end

%#ok<*AGROW,*HIST> This avoids code warnings for vectors and matrices extended in a loop
