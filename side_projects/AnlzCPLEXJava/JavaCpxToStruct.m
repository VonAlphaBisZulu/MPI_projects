function [scpx, varType] = JavaCpxToStruct(jcpx)
varType = class(jcpx);
switch varType
    %% CPLEX model
    case 'IloCplexExt'
        % Iterator is needed to extract information from model
        iter = jcpx.iterator;
        iR = 1;
        iT = 1;
        while iter.hasNext
            obj = iter.next;
            switch class(obj)
                case 'ilog.cplex.CpxLPMatrix' % The LP Matrix
                    vars    = obj.getNumVars;
                    ranges  = obj.getRanges;
                    numv    = length(vars);
                    numr    = length(ranges);
                    
                    % Variables
                    scpx.varName = cell.empty(numv,0);
                    scpx.varType = cell.empty(numv,0);
                    scpx.varLB   = double.empty(numv,0);
                    scpx.varUB   = double.empty(numv,0);
                    % Ranges and LP Matrix
                    scpx.ranName = cell.empty(numr,0);
                    scpx.ranType = cell.empty(numr,0);
                    scpx.ranLB   = double.empty(numr,0);
                    scpx.ranUB   = double.empty(numr,0);
                    scpx.LPmat   = zeros(numv,numr);

                    for i = 1:numv
                        [scpx.varName(i), scpx.varType(i), scpx.varLB(i), scpx.varUB(i)] = getVarNamesAndBounds(vars(i));
                    end
                    for i = 1:numr
                        [var, coeff] = getVarsAndCoeffFromExpr(ranges(i).getExpr);
                        for j = 1:length(var)
                            scpx.LPmat(strcmp(scpx.varName,var(j)),i) = coeff(j);
                        end
                        scpx.ranName(i) =  {char(ranges(i).getName)};
                        scpx.ranType(i) =  {char(ranges(i).getType)};
                        scpx.ranLB(i)   =  ranges(i).getLB;
                        scpx.ranUB(i)   =  ranges(i).getUB;
                    end
                case 'ilog.cplex.CpxRange'  % A supplementary constraint
                    if iR == 1
                        % Supplementary Constraints
                        scpx.conName = cell.empty(numr,0);
                        scpx.conType = cell.empty(numr,0);
                        scpx.conLB   = double.empty(numr,0);
                        scpx.conUB   = double.empty(numr,0);
                    end
                    scpx.conName(iR) = {char(obj.getName)};
                    scpx.conType(iR) = {char(obj.getType)};
                    scpx.conLB(iR)   = obj.getLB;
                    scpx.conUB(iR)   = obj.getUB;
                    
                    [var, coeff] = getVarsAndCoeffFromExpr(obj.getExpr);
                    con(iR,:) = [{var},coeff];
                    iR = iR+1;
                case 'ilog.cplex.CpxIfThen' % Logical and indicator constraints
                    if iT == 1
                        % If-Then conditions / indicator constraints
                        scpx.lgcName = cell.empty(numr,0);
                        scpx.lgcType = cell.empty(numr,0);
                        scpx.lgcLB   = double.empty(numr,0);
                        scpx.lgcUB   = double.empty(numr,0);
                    end
                    scpx.lgcName(iT) = {char(obj.getName)};
                    scpx.lgcType(iT) = {char(obj.getType)};
                    scpx.lgcLB(iT)   = obj.getLB;
                    scpx.lgcUB(iT)   = obj.getUB;
                    iT = iT+1;
                otherwise
                    error(['untreated class: ''' class(obj) '''']);
            end
        end       
        
    %% LP Matrix
    case 'ilog.cplex.CpxLPMatrix'
        vars    = jcpx.getNumVars;
        ranges  = jcpx.getRanges;
        numv    = length(vars);
        numr    = length(ranges);
        
        scpx.varName = cell.empty(numv,0);
        scpx.varType = cell.empty(numv,0);
        scpx.varLB   = double.empty(numv,0);
        scpx.varUB   = double.empty(numv,0);
        
        scpx.ranName = cell.empty(numr,0);
        scpx.ranType = cell.empty(numr,0);
        scpx.ranLB   = double.empty(numr,0);
        scpx.ranUB   = double.empty(numr,0);
        
        scpx.LPmat   = zeros(numv,numr);
        
        for i = 1:numv
            [scpx.varName(i), scpx.varType(i), scpx.varLB(i), scpx.varUB(i)] = getVarNamesAndBounds(vars(i));
        end
        for i = 1:numr
            [var, coeff] = getVarsAndCoeffFromExpr(char(ranges(i)));
            for j = 1:length(var)
                scpx.LPmat(strcmp(scpx.varName,var(j)),i) = coeff(j);
            end
            scpx.ranName(i) =  {char(ranges(i).getName)};
            scpx.ranType(i) =  {char(ranges(i).getType)};
            scpx.ranLB(i)   =  ranges(i).getLB;
            scpx.ranUB(i)   =  ranges(i).getUB;
        end
    %% equation or range
    case 'ilog.cplex.CpxRange'
        scpx.name = {char(jcpx.getName)};
        scpx.expr = char(jcpx.getExpr);
        [vars, scpx.coeff] = getVarsAndCoeffFromExpr(jcpx.getExpr);
        arrayfun(@(x) {char(x.getName)},vars);
    case {'ilog.cplex.CpxLinkedExpr', 'ilog.cplex.CpxQLNumExpr'}
        scpx.expr = char(jcpx);
        [vars, scpx.coeff] = getVarsAndCoeffFromExpr(jcpx);
        arrayfun(@(x) {char(x.getName)},vars);
    %% variable vector
    case {'ilog.concert.IloIntVar[]', 'ilog.concert.IloNumVar[]'}
        numv = length(jcpx);
        % init struct
        scpx.name = cell.empty(numv,0);
        scpx.varType = cell.empty(numv,0);
        scpx.LB = double.empty(numv,0);
        scpx.UB = double.empty(numv,0);
        % fill struct
        for i = 1:numv
            [scpx.name(i), scpx.varType(i), scpx.LB(i), scpx.UB(i)] = getVarNamesAndBounds(jcpx(i));
        end
    %% single variabe
    case 'ilog.cplex.CpxNumVar'
    %% unknown class
    otherwise
        error(['class ''' class(jcpx) ''' not covered by function']);
end
%% check if bounds are actually infinite
%  realmax==-scpx.ranLB
%  ...
        
function [name, vtype, lb, ub] = getVarNamesAndBounds(var)
    name    = {char(var.getName)};
    vtype   = {char(var.getType)};
    lb      = var.getLB;
    ub      = var.getUB;
end

function [var, coeff] = getVarsAndCoeffFromExpr(expr)
    linit = expr.linearIterator;
    qudit = expr.quadIterator;
    l = 1;
    while linit.hasNext
        var(l) = linit.next;
        coeff(l) = linit.getValue;
        l = l+1;
    end
    if qudit.hasNext
        error('this function does not cover quadratic expressions');
    end
end

function [] = interpretIfThen(IfThConst)
end
end
