classdef ElementaryVectorsEnumerator2 < CplexJava
    properties (SetAccess = 'protected')
        first_z
        err
        sol_var_start %A Java index of first solution variable in mipmat
        sol_var_num %A number of solution variables in mipmat
    end
    
    methods (Static)
        function max_dev= validate_z_vars2(zv)
            zv= unique(zv(:));
            if length(zv) > 1 && (zv(1) ~= 0 || zv(2) ~= 1)
                sel= round(zv) == 1;
                zv(sel)= zv(sel) - 1;
                max_dev= max(abs(zv));
                fprintf('Some Boolean variables assumed non-integer values, largest deviation was: %e\n', max_dev);
            end
        end
        function c = newline()
            c = char(10);
        end
    end
    
    methods (Access = 'public')
        function obj= ElementaryVectorsEnumerator2
            obj= obj@CplexJava();
            obj.err= false;
            obj.cpx.setOut([]);
        end
        
        function [obj, new_sols, finished]= calculate_further_solutions2(obj, enum_method, rem_evs, rem_time,inverseCutCount)
            if nargin < 5
                inverseCutCount = [];
            end
            newline = char(10);
            finished= false;
            cplex_inner= setup_cplex_inner_class_access();
            rem_time= max(rem_time, 1); %A continue for at least 1 second
            if ~isinf(rem_time)
                obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, rem_time);
            end
            
            switch enum_method
                case {1, 3} %A 1: solve to optimality; 3: return feasible solution if one exists after solve stops
                    try
                        fprintf('.');
                        obj.cpx.solve();
                        stat= obj.cpx.getStatus();
                        if obj.cpx.getCplexStatus().equals(cplex_inner.CplexStatus.AbortTimeLim)
                            if enum_method == 3 && stat.equals(cplex_inner.Status.Feasible)
                                new_sols= retrieve_optimal_solution2(obj);
                                fprintf('.');
                            else
                                new_sols= [];
                            end
                        else
                            if stat.equals(cplex_inner.Status.Optimal)
                                new_sols= retrieve_optimal_solution2(obj);
                                fprintf('.');
                            elseif stat.equals(cplex_inner.Status.Infeasible)
                                fprintf('\nNo further solutions found; minimal solution size is %d.\n', obj.ev_size_lb);
                                new_sols= [];
                                finished= true;
                            else
                                error('\nUnexpected status %s after solve()!\n', char(stat.toString()));
                            end
                        end
                    catch lasterror
                        lerr= lasterror;
                        fprintf('\nSevere error during solve(): %s\n', lerr.message);
                        new_sols= [];
                        obj.err= true;
                    end
                case 2 %A populate single cardinality
                    obj.evs_sz.setBounds(obj.ev_size_lb, obj.ev_size_lb);
                    if isinf(rem_evs)
                        obj.cpx.setParam(cplex_inner.IntParam.SolnPoolCapacity, 2.1e9);
                        obj.cpx.setParam(cplex_inner.IntParam.PopulateLim, 2.1e9);
                    else
                        obj.cpx.setParam(cplex_inner.IntParam.PopulateLim, rem_evs);
                    end
                    try
                        fprintf('Enumerating elementary vectors of size %d... \n', obj.ev_size_lb);
                        fprintf('%d Reactions are proposed to be added to model \n', length(inverseCutCount));
                        disp(['Populate..' newline]);
                        obj.cpx.populate();
                        limit_reached= obj.cpx.getCplexStatus().equals(cplex_inner.CplexStatus.AbortTimeLim)...
                            || obj.cpx.getCplexStatus().equals(cplex_inner.CplexStatus.PopulateSolLim);
                        stat= obj.cpx.getStatus();
                        if stat.equals(cplex_inner.Status.Optimal)
                            %A status optimal is reached as soon as the first solution is
                            %A found because there is no objective function during populate
                            [new_sols, num_sols]= retrieve_solutions_from_pool2(obj);
                            fprintf('Found %d of size %d.\n', num_sols, obj.ev_size_lb);
                            if ~limit_reached
                                obj.ev_size_lb= obj.ev_size_lb + 1;
                            end
                        elseif stat.equals(cplex_inner.Status.Infeasible)
                            new_sols= [];
                            if ~limit_reached
                                obj.ev_size_lb= obj.ev_size_lb + 1;
                            end
                            fprintf(' none found.\n');
                        elseif ~limit_reached
                            error('\nUnexpected status %s after populate()!\n', char(stat.toString()));
                        else %A limit_reached with status neither optimal nor infeasible
                            new_sols= [];
                        end
                    catch lasterror
                        lerr= lasterror;
                        fprintf('\nSevere error during solve(): %s\n', lerr.message);
                        obj.err= true;
                        new_sols= [];
                    end
                case 4
                    try
                        cplex_inner= setup_cplex_inner_class_access();  %P
                        obj.cpx.setParam(cplex_inner.IntParam.MIPSearch,2);
                        obj.cpx.solve();
                        stat= obj.cpx.getStatus();
                        if obj.cpx.getCplexStatus().equals(cplex_inner.CplexStatus.AbortTimeLim)
                            if enum_method == 4 && stat.equals(cplex_inner.Status.Feasible)
                                new_sols= retrieve_optimal_solution2(obj);
                                fprintf('.');
                            else
                                new_sols= [];
                            end
                        else
                            if stat.equals(cplex_inner.Status.Optimal)
                                new_sols= retrieve_optimal_solution2(obj);
                                fprintf('.');
                            elseif stat.equals(cplex_inner.Status.Infeasible)
                                fprintf('\nNo further solutions found; minimal solution size is %d.\n', obj.ev_size_lb);
                                new_sols= [];
                                finished= true;
                            else
                                error('\nUnexpected status %s after solve()!\n', char(stat.toString()));
                            end
                        end
                    catch lasterror
                        lerr= lasterror;
                        fprintf('\nSevere error during solve(): %s\n', lerr.message);
                        obj.err= true;
                        new_sols= [];
                    end
            end % switch
        end
        
        function obj= add_evs_exclusion_constraints2(obj, zv, support_size)
            if support_size == 0
                disp('Empty exclusion constraint not added.');
            else
                m= obj.mipmat.getNrows();
                obj.mipmat.addRows(zeros(1, size(zv, 2)), support_size - 1, [], []);
                [i, j, val]= find(zv);
                obj.mipmat.setNZs(m + j - 1, obj.first_z + i - 1, val); %A Java indices, implicit transposition
            end
        end
        
        function new_sol= retrieve_optimal_solution2(obj)
            new_sol= obj.cpx.getValues(obj.mipmat, obj.sol_var_start, obj.sol_var_num);
        end
        
        function [new_sols, num_sols]= retrieve_solutions_from_pool2(obj)
            num_sols= obj.cpx.getSolnPoolNsolns();
            new_sols= zeros(obj.sol_var_num, num_sols);
            for i= 1:num_sols
                new_sols(:, i)= obj.cpx.getValues(obj.mipmat, obj.sol_var_start, obj.sol_var_num, i-1);
            end
        end
        
        function obj= set_objective_function2(obj)
            % only sets the objective function if none is present or if the
            % current one is a constant
            if isempty(obj.cpx.getObjective())
                obj.cpx.addMinimize().setExpr(obj.cpx.sum([obj.z_vars obj.iy_vars]));
            else
                objective= obj.cpx.getObjective();
                if ~objective.getExpr().linearIterator().hasNext()
                    objective.setExpr(obj.cpx.sum([obj.z_vars obj.iy_vars]));
                end
            end
        end
        
        function [obj, obj_expr]= clear_objective_function2(obj)
            objective= obj.cpx.getObjective();
            obj_expr= objective.getExpr();
            objective.clearExpr();
        end
        
        function obj= setup_fixed_cplex_parameters2(obj, enum_method)
            cplex_inner= setup_cplex_inner_class_access();
            try
                numcores = feature('numcores');
                obj.cpx.setParam(cplex_inner.IntParam.Threads, numcores);
            catch lasterror
                warning(['could not find number of cores, uses CPLEX standard', lasterror.message])
            end
            obj.cpx.setParam(cplex_inner.IntParam.ParallelMode, 1); %
            if ispc
               ret = getenv('COMPUTERNAME');
            else      
               ret = getenv('HOSTNAME');     
            end
            if isempty(ret)
                [st,ret] = system('uname -n');
            end
            switch enum_method % set up method-specific parameters that do not change during enumeration
                case {1, 3, 4} %A single solve
                    %A setting emphasis to feasibility appears to be faster than optimality
                    obj.cpx.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Feasibility);
                    obj.cpx.setParam(cplex_inner.IntParam.MIPDisplay, 0);
                    %P Barrier mode causes crashes, so dual simplex mode is
                    % chosen as fastest alternative
                    % if running on node: 
                    if ~isempty(strfind(ret,'node')) && ~isempty(strfind(ret,'perceus'))
                        obj.cpx.setParam(cplex_inner.IntParam.RootAlg, 0);
                    else
                        obj.cpx.setParam(cplex_inner.IntParam.RootAlg, 2);
                    end
                case 2 %A populate single cardinality
                    obj.cpx.setParam(cplex_inner.IntParam.SolnPoolIntensity, 4);
                    obj.cpx.setParam(cplex_inner.IntParam.MIPDisplay, 0);
                    obj.cpx.setParam(cplex_inner.IntParam.MIPEmphasis, cplex_inner.MIPEmphasis.Feasibility);
                    obj.cpx.setParam(cplex_inner.DoubleParam.SolnPoolAGap, 0);
                    %A deterministic parallel mode without objective function using branch & cut appears to be most efficient
                    obj.cpx.setParam(cplex_inner.IntParam.ParallelMode, 1);
                    obj.cpx.setParam(cplex_inner.IntParam.MIPSearch, 1);
                    %P Barrier mode causes crashes, so dual simplex mode is
                    % chosen as fastest alternative
                    if ~isempty(strfind(ret,'node')) && ~isempty(strfind(ret,'perceus'))
                        obj.cpx.setParam(cplex_inner.IntParam.RootAlg, 0);
                    else
                        obj.cpx.setParam(cplex_inner.IntParam.RootAlg, 2);
                    end
            end
        end
    end % methods
end
