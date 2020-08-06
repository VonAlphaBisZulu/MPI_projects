function [obj, mcs, status, comp_time]= shortest_minimal_cut_sets2(obj, max_ev_size, num_evs_lim, time_limit, enum_method, validate, cplex_parset,inverseCutCount)
% enum_method:
% 1: solve
% 2: populate
% 3: solve to see if a feasible solution exists
% 4: solve to see if a feasible solution exists and return n solutions

if nargin < 6 
  validate= false;
  if nargin < 5
    enum_method= 2;
    if nargin < 4 || isempty(time_limit)
      time_limit= Inf;      
    end
  end
else
    if isempty(validate)
        validate= false;
    end
end
if nargin < 8 
    inverseCutCount = [];
end

n= obj.sol_var_num/2;
mcs= false(n, 0);
rem_time= Inf;
start_time= clock();
status=[];
newline = char(10);

enum_method= fix(enum_method); % make integer
if enum_method < 1 || enum_method > 4
  fprintf('Invalid enumeration method; exiting.');
  obj.err= true;
  return
end
if enum_method == 3
  num_evs_lim= 1;
mcs= false(n, 0);
  max_ev_size= n;
end

if nargin < 7
  obj= setup_fixed_cplex_parameters2(obj, enum_method);
else
    if ~isempty(cplex_parset)
        obj.cpx.setParameterSet(cplex_parset);
    else
        obj= setup_fixed_cplex_parameters2(obj, enum_method);
    end
end
if enum_method == 1 %A solve
  obj= set_objective_function2(obj);
end
num_new_evs= 0; 
continue_loop= true;
disp('Starting enumeration of minimal cut sets.');
disp(['Solve..' newline]);
while continue_loop && obj.ev_size_lb <= max_ev_size
  rem_evs= num_evs_lim - num_new_evs;
  if rem_evs <= 0
    disp('Solution limit reached.');
    break;
  end
  if ~isinf(time_limit)
    rem_time= time_limit - etime(clock(), start_time);
    if rem_time <= 0
      disp('Time limit reached.');
      break;
    end
  end
   [obj, new_sols, finished]= calculate_further_solutions2(obj, enum_method, rem_evs, rem_time,inverseCutCount);
  status= obj.cpx.getStatus(); % remember the status before modifying the MILP
  if finished || obj.err
    continue_loop= false;
  end
  if ~isempty(new_sols)
      try
          ElementaryVectorsEnumerator2.validate_z_vars2(new_sols);
          new_sols= logical(round(new_sols));
          inverseCountCuts = ones(2*n,size(new_sols,2));
          inverseCountCuts([inverseCutCount;n+inverseCutCount],:) = -1 ;
          support_size= sum(new_sols.*inverseCountCuts, 1);
          real_support_size = sum(new_sols);
          if enum_method ~= 1 || support_size <= max_ev_size
              new_sols= new_sols(1:n, :) | new_sols(n+1:end, :);
              mcs= [mcs, new_sols];
              num_new_evs= num_new_evs + size(new_sols, 2);
              if enum_method ~= 3 %A only integrate exclusion constraints for optimal solutions
                  obj= add_mcs_exclusion_constraints2(obj, new_sols, real_support_size);
              end
          end
          if enum_method ~= 3 && enum_method ~= 4 && support_size(1) > obj.ev_size_lb
              obj.ev_size_lb= support_size(1);
              fprintf('\nIncreasing minimal set size to %d\n', obj.ev_size_lb);
              obj.evs_sz.setLB(obj.ev_size_lb);
          end
      catch warn
        warning(warn.message);
        disp('the mcs found up to this point will be returned.');
        continue_loop = false;
      end
  end
end % while
comp_time= etime(clock(), start_time);
fprintf('\nFinished\n');
