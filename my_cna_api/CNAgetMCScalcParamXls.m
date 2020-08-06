function [T, t, D, d,notknockable,reacMin,reacMax,xors] = CNAgetMCScalcParamXls( cnap, xls_input_filename)

% bIgnFE (ignore false entries): If entries cannot be found, a warning is thrown
% instead of an error - useful if genome scale setups are used to calculate
% cut sets in a core model

reacMin = cnap.reacMin;
reacMax = cnap.reacMax;

%% read out cMCS-Calculation parameters
TableReads = loadSpecReacXLStoStrArray(xls_input_filename);

% ------------------ identify constraint sheet ---------------------------
% Find "sheet" with parameters for cMCS calculation
[nkorow,nkocol,cMCSSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'notknockable')));

[tRrow,tRcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'targetR'));
[dRrow,dRcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'desiredR'));
[gCrow,gCcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'generalConst'));
[RMinrow,RMinCcol,~] = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmin'));
[RMaxrow,RMaxcol,~]  = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmax'));
[kirow,kicol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'knockin'));
[xorrow,xorcol,~]    = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'xor'));

% ================ Target Reactions ===================
for i = 1:length(tRrow)
    lastrow = tRrow(i)+find(strcmp(strtrim(TableReads(tRrow(i):end,tRcol(i),cMCSSheet)),''),1,'first')-2;
    if isempty(lastrow)
        lastrow = tRrowtRrow(i)+find(~strcmp(strtrim(TableReads(tRrowtRrow(i):end,tRcol(i),cMCSSheet)),''),1,'last')-1;
    end
    targetReacs{i} = [cellstr(TableReads(tRrow(i)+1:lastrow,tRcol(i):tRcol(i)+1,cMCSSheet)) num2cell(str2double(TableReads(tRrow(i)+1:lastrow,tRcol(i)+2,cMCSSheet)))];
end

% ================ Desired Reactions ===================
for i = 1:length(dRrow)
    lastrow = dRrow(i)+find(strcmp(strtrim(TableReads(dRrow(i):end,dRcol(i),cMCSSheet)),''),1,'first')-2;
    if isempty(lastrow)
        lastrow = dRrow(i)+find(~strcmp(strtrim(TableReads(dRrow(i):end,dRcol(i),cMCSSheet)),''),1,'last')-1;
    end
    desiredReacs{i} = [cellstr(TableReads(dRrow(i)+1:lastrow,dRcol(i):dRcol(i)+1,cMCSSheet)) num2cell(str2double(TableReads(dRrow(i)+1:lastrow,dRcol(i)+2,cMCSSheet)))];
end

% ================ general Constrains ===================
lastrow = gCrow+find(strcmp(strtrim(TableReads(gCrow:end,gCcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = gCrow+find(~strcmp(strtrim(TableReads(gCrow:end,gCcol,cMCSSheet)),''),1,'last')-1;
end
genericConstrains = [cellstr(TableReads(gCrow+1:lastrow,gCcol:gCcol+1,cMCSSheet)) num2cell(str2double(TableReads(gCrow+1:lastrow,gCcol+2,cMCSSheet)))];

targetReacs  = cellfun(@(x) [x;  genericConstrains],targetReacs,'UniformOutput' , false);
desiredReacs = cellfun(@(x) [x;  genericConstrains],desiredReacs,'UniformOutput' , false);

% ================ Rmin ===================
lastrow = RMinrow+find(strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMinrow+find(~strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'last')-1;
end
Rmin = [cellstr(TableReads(RMinrow+1:lastrow,RMinCcol,cMCSSheet)) num2cell(str2double(TableReads(RMinrow+1:lastrow,RMinCcol+1,cMCSSheet)))];

% ================ Rmax ===================
lastrow = RMaxrow+find(strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMaxrow+find(~strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'last')-1;
end
Rmax = [cellstr(TableReads(RMaxrow+1:lastrow,RMaxcol,cMCSSheet)) num2cell(str2double(TableReads(RMaxrow+1:lastrow,RMaxcol+1,cMCSSheet)))];

% ================ notknockable ===================
lastrow = nkorow+find(strcmp(strtrim(TableReads(nkorow:end,nkocol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = nkorow+find(~strcmp(strtrim(TableReads(nkorow:end,nkocol,cMCSSheet)),''),1,'last')-1;
end
nokoable = cellstr(TableReads(nkorow+1:lastrow,nkocol,cMCSSheet));
notknockable = [];
% ================ knockin ===================
lastrow = kirow+find(strcmp(strtrim(TableReads(kirow:end,kicol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = kirow+find(~strcmp(strtrim(TableReads(kirow:end,kicol,cMCSSheet)),''),1,'last')-1;
end
kin = cellstr(TableReads(kirow+1:lastrow,kicol,cMCSSheet));
knockinable = [];
% ================ xor ===================
for i = 1:length(xorrow)
    lastrow = xorrow(i)+find(strcmp(strtrim(TableReads(xorrow(i):end,xorcol(i),cMCSSheet)),''),1,'first')-2;
    if isempty(lastrow)
        lastrow = xorrow(i)+find(~strcmp(strtrim(TableReads(xorrow(i):end,xorcol(i),cMCSSheet)),''),1,'last')-1;
    end
    xor(i) = {cellstr(TableReads(xorrow(i)+1:lastrow,xorcol(i),cMCSSheet))};
end

%% identify not knockable reactions

for ink = 1:length(nokoable)
    if any(strcmp(strtrim(cellstr(cnap.reacID)), nokoable(ink)))
        notknockable = [notknockable find(strcmp(strtrim(cellstr(cnap.reacID)), nokoable(ink)))];
    else
        warning(['Reaction ' char(nokoable(ink)) ' could not be found and wasn''t set to ''not knockable''.']);
    end
end

%% identify reactions that can be knocked in

for ink = 1:length(kin)
    if any(strcmp(strtrim(cellstr(cnap.reacID)), kin(ink)))
        knockinable = [knockinable find(strcmp(strtrim(cellstr(cnap.reacID)), kin(ink)))];
    else
        warning(['Reaction ' char(kin(ink)) ' could not be found and isn''t considered as ''knock in candidate''.']);
    end
end

% if knockins are contained in 'notknockable' these reactions are removed from there
% because reactions must be knockable - even to be activated
notknockable(~ismember(notknockable,knockinable));

%% identify reactions that exclude each other
if ~isempty(xorrow)
    xors = repmat({double.empty(1,0)},length(xor),1);
    for x = 1:length(xor)
        for i = 1:length(xor{x})
            if any(strcmp(strtrim(cellstr(cnap.reacID)), xor{x}(i)))
                xors{x} = [xors{x} find(strcmp(strtrim(cellstr(cnap.reacID)), xor{x}(i)))];
            else
                warning(['Reaction ' char(kin(ink)) ' could not be found and isn''t considered as ''knock in candidate''.']);
            end
        end
    end
end


%% set Rmin and Rmax

for irmn = 1:length(Rmin)
    if any(strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))
        reacMin((strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))) = cell2mat(Rmin(irmn,2));
    else
        warning(['Reaction ' char(Rmin(irmn,1)) ' could not me found. Boundaries were not set.']);
    end
end

for irmx = 1:length(Rmax)
    if any(strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))
        reacMax((strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))) = cell2mat(Rmax(irmx,2));
    else
        warning(['Reaction ' char(Rmax(irmx,1)) ' could not me found. Boundaries were not set.']);
    end
end

%% generate t, T, d, D
% generate t, T
for i = 1:length(targetReacs)
    [T{i},t{i}] = genV(targetReacs{i},cellstr(cnap.reacID),reacMin,reacMax);
end

% generate d, D
for i = 1:length(desiredReacs)
    [D{i},d{i}] = genV(desiredReacs{i},cellstr(cnap.reacID),reacMin,reacMax);
end

function [V,v] = genV(constraints,reacID,rMin,rMax)
% Generate Vectors V and v so that V*r <= v
% input: some constraint seperated 3 three cells. e.g.:
%           r_1 + r_4 / r_3 - r_2    |    >=    |    a

    V = [];
    v = [];

    for j = 1:size(constraints,1)
        % get right hand side
        a = constraints{j,3};
        % get direction of inequality
        switch constraints{j,2}
            case '<='
                eqop = 1;
            case '>='
                eqop = -1;
            otherwise
                error('please define inequality either as ''<='' or ''>=''');
        end
        % split into numerator an divisor and get coefficients
        numDiv = split(constraints{j,1},'/');
        % get variables and coefficients for vector
        [num, cnum] = findReacAndCoeff(numDiv(1),reacID);
            % fractional constraint ispreprocessed
        if length(numDiv) == 2
            [div, cdiv] = findReacAndCoeff(numDiv(2),reacID);
            % check if all reactions take identical signs (+ or -)
            % this is needed to do the equation rearrangement. If the signs
            % are ambigous, a lot of case differentiations would be needed,
            % what is not done here.
            if any(rMin(div).*rMax(div) < 0)
                error(['reactions that are part of the divisor inside the '... 
                        'target constraints must not span a positive AND negative range']);
            else
                rDir = sign(rMin(div)+rMax(div));
                if any(rDir.*cdiv' > 0) &&  any(rDir.*cdiv' < 0)
                    error(['reactions that are part of the divisor inside the '...
                           'target constraints must all have the same direction '...
                           '(positive or negative). Please check if coefficients and '...
                           'reaction ranges lead to all positive or all negative variables.']);
                else
                    if sign(sum(rDir.*cdiv')) == -1 % if the divisor is all negative
                        % change direction of inequality
                        eqop = -eqop;
                    end
                end
            end
            cdiv = -a*cdiv; % transformation to following form:
            a = 0;
            % (cnum num) - a*cdiv1*div1 - a*cdiv2*div2 <= 0
            num  = [num,   div];
            cnum = [cnum, cdiv];
        end
        % constraint is generated
        if eqop == -1  % if inequality is >= the whole term is multiplied with -1
            cnum = -cnum;
            a    = -a;
        end
        v(j,1)   = a;
        V(j,:) = full(sparse(num,1,cnum,length(reacID),1));
    end

    function [ridx,coeff] = findReacAndCoeff(eq,reacID)
        % extract reac indices and coefficients from equation
        r = regexp(eq, reacID, 'match');
        ridx = find(~cellfun(@isempty,r))';
        r = [r{ridx}];
        for k = 1:length(r)
            c = regexp(eq, ['(\s|\d|-|\.)*?(?=' r{k} ')'], 'match');
            c = strrep(char(c{:}),' ','');
            switch c
                case ''
                    coeff(k) = 1;
                case '-'
                    coeff(k) = -1;
                otherwise
                    coeff(k) = str2num(c);
            end
        end
    end
end
end