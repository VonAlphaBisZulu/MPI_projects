classdef ConstrainedMinimalCutSetsEnumerator2 < ElementaryVectorsEnumerator2
  properties (SetAccess = 'immutable')
    flux_vars
    evs_sz %
    z_vars
    iy_vars
  end
  
  properties
    ev_size_lb
  end

%% Constructor
  methods (Access = 'public')                                                          % | params for integrated mode | if some cuts are count negatively 
    function obj= ConstrainedMinimalCutSetsEnumerator2(st, irr, targets, inh, ubi, cuts,   flux_lb, flux_ub, des, db, inverseCutCount)
      if nargin < 7 % MISSING: flux_lb, flux_ub, des, db, inverseCutCount
        include_mnet= false;
        if nargin < 6
          cuts= [];
          if nargin < 5
            inh= [];
            ubi= [];
          end
        end
      elseif nargin >= 7
        if ~isempty(flux_lb) && ~isempty(flux_lb)
            include_mnet= true; % integrated mode with flux bounds for all
                                % reactions that are KOable:
                                % -Inf <= R + ZP * ub <= ub
                                % -Inf <= R + ZN * ub <= ub
                                %   lb <= R + ZP * lb <= Inf
                                %   lb <= R + ZN * lb <= Inf
        else
            include_mnet= false;
        end
      end
      if nargin < 11
          inverseCutCount = []; % list with inversely count reactions
      end
      
      [m, n]= size(st);
      st= sparse(st);
      irr= logical(irr);
      if ~isempty(targets) && ~isempty(ubi)
        error('Choose either target reactions or inhomogeneous constraints.');
      end
      if ~isempty(targets) %A# set up targets column-wise
        if isvector(targets)
          if any(targets & ~irr)
            error('Can only select irreversible reactions as targets');
          end
          targets= targets(:);
        else
          if any(any(targets, 1) & ~irr)
            error('Can only select irreversible reactions as targets');
          end
          targets= targets';
        end
      else
        targets= [];
      end
      if isempty(inh)
        inh= [];
      end
      if ~isempty(ubi)
        ubi= ubi(:)';
      end
      ic= size(inh, 1);
      tc= size(targets, 2); %A# T was transposed above
      numsn= nums(n);
      obj= obj@ElementaryVectorsEnumerator2();
      % variable names, except z, according to the bioinformatics paper:
      % v (2*n) >= 0, q (1) > 0, w (ic + tc) >= 0, u (m), z (2*n) [0,1]
      nv= 2*n + 1 + ic + tc + m;
      vqwu_vars= obj.cpx.numVarArray(nv, [zeros(1, 2*n), 1, zeros(1, ic + tc), -Inf(1, m)], Inf(nv, 1),...
        [strcat('VP', numsn); strcat('VN', numsn); 'Q';strcat('W', nums(ic+tc)); strcat('U', nums(m))]);
      obj.mipmat.addCols(vqwu_vars);
      if include_mnet
        obj.flux_vars= obj.cpx.numVarArray(n, flux_lb, flux_ub, strcat('R', numsn));
        obj.mipmat.addCols(obj.flux_vars);
        first_fv= int32(nv); % Java index
        obj.first_z= nv + n; % Java index
      else
        obj.flux_vars= [];
        obj.first_z= nv; % Java index
      end
      obj.sol_var_start= obj.mipmat.getNcols();
      obj.sol_var_num= 2*n;
      
      nonInv = cellstr(num2str(setdiff(str2num(char(numsn)),inverseCutCount)));
      if ~isempty(inverseCutCount)
        inv    = cellstr(num2str(                             inverseCutCount));
      else
          inv = {};
      end
      ni  = length(inv);
      
      obj.z_vars= obj.cpx.boolVarArray(2*n, [strcat('ZP', numsn); strcat('ZN', numsn)]); %A# 1..n: z_pos, n+1..2*n: z_neg
      obj.mipmat.addCols(obj.z_vars); %A put z_vars at the end of mipmat
      for i= n + find(irr) %A quick fix: disable z_neg of irreversible reactions
        obj.z_vars(i).setUB(0);
      end
      if ~isempty(cuts)
        cuts= cuts(:)'; %A make row vector
        for i= find(~cuts)
          obj.z_vars(i).setUB(0);
          obj.z_vars(i+n).setUB(0);
        end
      end
      
      obj.iy_vars= obj.cpx.boolVarArray(ni, strcat('IY', inv));
      if ni
        obj.mipmat.addCols(obj.iy_vars);
      end
      
      %P make z_var-List that doesn't contain the inverse cuts
      obj_expr= obj.cpx.sum([obj.iy_vars,obj.z_vars([str2num(char(nonInv)) n+str2num(char(nonInv))])]);
      %obj_expr= obj.cpx.sum(obj.z_vars); %A# MCS size
      obj.evs_sz= obj.cpx.range(0, obj_expr, n); %A# this constraint controls the MCS size
      obj.mipmat.addRow(obj.evs_sz); %A# put this in the first row of mipmat
%       obj.cpx.addMinimize().setExpr(obj_expr);
      obj.ev_size_lb= 0;
      
      ub= zeros(n, 1);
      ub(irr)= Inf;
      obj.mipmat.addRows(zeros(n, 1), ub, [], []);
      [i, j, val]= find([inh', -targets, st']);
      % i instead (i-1) because the first row of mipmat is mcs_sz; j instead of
      % (j-1) because of the q column:
      obj.mipmat.setNZs(i, j + 2*n, val);
      obj.mipmat.setNZs(1:n, 0:(n-1), ones(1, n)); % 1:n instead of 0:(n-1) for the same reason
      obj.mipmat.setNZs(1:n, n:(2*n-1), -ones(1, n));
      
      ind= find(ubi);
      val= [1, ubi(ind), -ones(1, tc)];
      ind= int32([1, ind + 1, length(ind)+2:length(ind)+1+tc] + 2*n - 1); %A# Java indices
      obj.mipmat.addRow(0, 0, ind, val);
      
      if include_mnet
        obj.mipmat.addRows([zeros(m, 1); -Inf(length(db), 1)], [zeros(m, 1); db], [], []);
        [i, j, val]= find([st; des]);
        obj.mipmat.setNZs(i + n + 1, first_fv + int32(j - 1), val);
      end
      indic= javaArray('ilog.concert.IloConstraint', 4*n);
      %A setting up indicators using ranges from LPMatrices is faster
      %A than using expressions over variables
      zeq0= obj.cpx.LPMatrix();
      zeq0.addCols(obj.z_vars);
      zeq0.addRows(zeros(2*n, 1), zeros(2*n, 1), [], []);
      zeq0.setNZs(0:(2*n)-1, 0:(2*n)-1, ones(1, 2*n));
      veq0= obj.cpx.LPMatrix();
      veq0.addCols(vqwu_vars(1:2*n));
      veq0.addRows(zeros(2*n, 1), zeros(2*n, 1), [], []);
      veq0.setNZs(0:(2*n)-1, 0:(2*n)-1, ones(1, 2*n));      
      zeq1= obj.cpx.LPMatrix();
      zeq1.addCols(obj.z_vars);
      zeq1.addRows(ones(2*n, 1), ones(2*n, 1), [], []);
      zeq1.setNZs(0:(2*n)-1, 0:(2*n)-1, ones(1, 2*n));
      vge1= obj.cpx.LPMatrix();
      vge1.addCols(vqwu_vars(1:2*n));
      vge1.addRows(ones(2*n, 1), Inf(2*n, 1), [], []);
      vge1.setNZs(0:(2*n)-1, 0:(2*n)-1, ones(1, 2*n));

      for i= 1:2*n
        indic(2*i-1)= obj.cpx.add(obj.cpx.ifThen(zeq0.getRange(i-1), veq0.getRange(i-1)));
        indic(2*i)= obj.cpx.add(obj.cpx.ifThen(zeq1.getRange(i-1), vge1.getRange(i-1)));
      end
      i= obj.mipmat.getNrows();

      obj.mipmat.addRows(zeros(n, 1), ones(n, 1), [], []); %A# z_pos(i) and z_neg(i) must not both be 1
      obj.mipmat.setNZs(i:i+n-1, obj.first_z:obj.first_z+n-1, ones(1, n));
      obj.mipmat.setNZs(i:i+n-1, obj.first_z+n:obj.first_z+2*n-1, ones(1, n));
      
      % inversely count cuts
      if ni
          yeq0= obj.cpx.LPMatrix();
          yeq0.addCols(obj.iy_vars);
          yeq0.addRows(zeros(ni, 1), zeros(ni, 1), [], []);
          yeq0.setNZs(0:(ni)-1, 0:(ni)-1, ones(1, ni));
          zeq0= obj.cpx.LPMatrix();
          zeq0.addCols(obj.z_vars(1:2*n));
          zeq0.addRows(zeros(ni, 1), zeros(ni, 1), [], []);
          zeq0.setNZs(repmat(0:(ni)-1,1,2), [str2num(char(inv))-1;str2num(char(inv))+n-1]', ones(1, 2*ni));      
          yeq1= obj.cpx.LPMatrix();
          yeq1.addCols(obj.iy_vars);
          yeq1.addRows(ones(ni, 1), ones(ni, 1), [], []);
          yeq1.setNZs(0:(ni)-1, 0:(ni)-1, ones(1, ni));
          zeq1= obj.cpx.LPMatrix();
          zeq1.addCols(obj.z_vars(1:2*n));
          zeq1.addRows(ones(ni, 1), ones(ni, 1), [], []);
          zeq1.setNZs(repmat(0:(ni)-1,1,2), [str2num(char(inv))-1;str2num(char(inv))+n-1]', ones(1, 2*ni));

          for i = 1:ni
            indic(2*i-1)= obj.cpx.add(obj.cpx.ifThen(yeq0.getRange(i-1), zeq1.getRange(i-1)));
            indic(2*i)= obj.cpx.add(obj.cpx.ifThen(yeq1.getRange(i-1), zeq0.getRange(i-1)));
          end
      end
      
      if include_mnet
        if isempty(cuts)
          ind= 1:n;
        else
          % for a currently unlcear reason it appears to be important to not link
          % the (fixed) z_vars of reactions which cannot be cut to their
          % corresponding flux_vars
          % kann das daran liegen daß bei essentiellen Reaktionen die bounds die
          % 0 ausschliessen und durch Hilfsvariablen für das ifThen
          % Inkonsistenzen entstehen?
          ind= find(cuts);
        end
%         tic;
        for i= ind
          if flux_ub(i) ~= 0
            %       obj.cpx.addLe(obj.flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.cpx.sum(obj.z_vars(i), obj.z_vars(i+n))), flux_ub(i)));
            % unclear whether a single or separate constraints lead to better performance
            obj.cpx.addLe(obj.flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.z_vars(i)), flux_ub(i)));
            obj.cpx.addLe(obj.flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.z_vars(i+n)), flux_ub(i)));
          end
          if flux_lb(i) ~= 0
            %       obj.cpx.addGe(flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.cpx.sum(obj.z_vars(i), obj.z_vars(i+n))), flux_lb(i)));
            % unclear whether a single or separate constraints lead to better performance
            obj.cpx.addGe(obj.flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.z_vars(i)), flux_lb(i)));
            obj.cpx.addGe(obj.flux_vars(i), obj.cpx.prod(obj.cpx.diff(1, obj.z_vars(i+n)), flux_lb(i)));
          end
        end
%         toc;
      end
    end
  
    function obj= add_mcs_exclusion_constraints2(obj, zv, support_size)
      if support_size == 0
        disp('Empty exclusion constraint not added.');
      else
        zv= [zv; zv];
        m= obj.mipmat.getNrows();
        obj.mipmat.addRows(zeros(1, size(zv, 2)), support_size - 1, [], []);
        [i, j, val]= find(zv);
        obj.mipmat.setNZs(m + j - 1, obj.first_z + i - 1, val); %A Java indices, implicit transposition
      end
    end
    
    % method declarations
    [obj, mcs, status, comp_time]= shortest_minimal_cut_sets2(obj, max_ev_size, num_evs_lim, time_limit, enum_method, validate, cplex_parset,inverseCutCount)
  end
end

function res= nums(n)
res= strtrim(cellstr(num2str((1:n)')));
end
