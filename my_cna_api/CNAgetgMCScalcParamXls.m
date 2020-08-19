function [T, t, D, d,rkoCost,rkiCost,reacMin,reacMax,gkoCost,gkiCost,idx] = CNAgetgMCScalcParamXls( cnap, xls_prod_filename, xls_subs_filename, genes)
if ~exist('genes','var')
    genes = [];
    gkoCost = [];
    gkiCost = [];
end

reacMin = cnap.reacMin;
reacMax = cnap.reacMax;

%% read out cMCS-Calculation parameters
table_prod = loadSpecReacXLStoStrArray(char(xls_prod_filename));

Rmin = readOut(table_prod,'Rmin',2);
Rmax = readOut(table_prod,'Rmax',2);
rkoCost_cell = readOut(table_prod,'ko_r',2);
rkiCost_cell = readOut(table_prod,'ki_r',2);

% reaction KO/KI costs
rkoCost = nan(1,cnap.numr);
for i = 1:size(rkoCost_cell,1)
    rkoCost(contains(cellstr(cnap.reacID),rkoCost_cell{i,1})) = str2double(rkoCost_cell{i,2});
end
rkiCost = nan(1,cnap.numr);
for i = 1:size(rkiCost_cell,1)
    rkiCost(contains(cellstr(cnap.reacID),rkiCost_cell{i,1})) = str2double(rkiCost_cell{i,2});
end

if ~isempty(genes)
    try
        gkoCost_cell = readOut(table_prod,'ko_g',2);
        gkiCost_cell = readOut(table_prod,'ki_g',2);
    catch
        gkoCost_cell = cell.empty(0,2);
        gkiCost_cell = cell.empty(0,2);
    end
    % gene KO/KI costs
    gkoCost = ones(1,length(genes)); % by default, all genes are knockable
    gkiCost = nan(1,length(genes)); %                 and not knock-inable
    if ((~isempty(gkoCost_cell) || ~isempty(gkoCost_cell)) && ~exist('genes','var'))
        warning('The information for gene knock-outs/ins was incomplete and gene-KOs will be ignored. Check if you provided a list of all gene names');
    else
        for i = 1:size(gkoCost_cell,1)
            gkoCost(contains(genes,gkoCost_cell{i,1})) = str2double(gkoCost_cell{i,2});
        end
        for i = 1:size(gkiCost_cell,1)
            gkiCost(contains(genes,gkiCost_cell{i,1})) = str2double(gkiCost_cell{i,2});
        end
    end
end

%% set Rmin and Rmax
Rmin_idx = findStrPos(cnap.reacID,Rmin(:,1));
reacMin(Rmin_idx(Rmin_idx~=0)) = str2double(Rmin(Rmin_idx~=0,2));
if any(Rmin_idx == 0)
    warning(['Reaction(s) ' strjoin(Rmin(Rmin_idx == 0,1),', ') ' could not be found. Boundaries were not set.']);
end
Rmax_idx = findStrPos(cnap.reacID,Rmax(:,1));
reacMax(Rmax_idx(Rmax_idx~=0)) = str2double(Rmax(Rmax_idx~=0,2));
if any(Rmax_idx == 0)
    warning(['Reaction(s) ' strjoin(Rmax(Rmax_idx == 0,1),', ') ' could not be found. Boundaries were not set.']);
end
cnap.reacMin = reacMin;
cnap.reacMax = reacMax;

for i = 1:length(xls_subs_filename)
    table_subst{i} = loadSpecReacXLStoStrArray(char(xls_subs_filename(i)));
end

%% generate t, T, d, D
% generate t, T
targetReacs  = readOut(table_prod,'targetR',3);
if exist('targetReacs','var') && ~iscell(targetReacs{1})
    targetReacs = {targetReacs};
end
desiredReacs = readOut(table_prod,'desiredR',3);
if exist('desiredReacs','var') && ~iscell(desiredReacs{1})
    desiredReacs = {desiredReacs};
end
for i = 1:length(xls_subs_filename)
    uptake_reac_name = readOut(table_subst{i},'substrate',1);
    yield_factor     = readOut(table_subst{i},'yield_factor',1);
    subst(i).fac_and_name = [yield_factor{:},' ',uptake_reac_name{:}];
    export_reac = readOut(table_subst{i},'active_export_reac',1);
    if ~isempty(export_reac)
        subst(i).subst_export_constr = [export_reac {'<='} {'0'}];
    else
        subst(i).subst_export_constr = cell.empty(0,3);
    end
    kicost  = str2double(readOut(table_subst{i},'kiCost',1));        % set KI-costs
    rkiCost(contains(cellstr(cnap.reacID),uptake_reac_name)) = kicost;
end
target_w_subst_placeholder = find(cellfun(@(x) any(any(contains(x,'#subst#'))),targetReacs));
% if no target region definition contains the keyowrd #subst#, even though substrate files were added, throw warning
if isempty(target_w_subst_placeholder) && isempty(xls_subs_filename)
    warning(['None of the target region definitions contains ''#subst#''',...
             'so substrates from separate excel files are ignored for the definition of target regions']);
end
% replace #subst# in target regions
for i = target_w_subst_placeholder % if one region definition contains #subst#, replace it / generate target regions
    comb = cellfun(@(x) double(isempty(x)):1, {subst(:).subst_export_constr},'UniformOutput',0);
    comb = combvec(comb{:}); % vector that determines which reactions are active in which target region
    comb = comb(:,any(comb,1)); % remove empty substrate sets
    new_target = repmat(targetReacs(i),1,size(comb,2));
    for j = 1:size(comb,2)
        subst_idx = find(comb(:,j));
        new_target{j} = [new_target{j}; vertcat(subst(subst_idx).subst_export_constr)]; % add export constraints
        subst_str = ['(' strjoin({subst(subst_idx).fac_and_name},' ') ')']; % build string for substrate combination
        new_target{j} = strrep(new_target{j},'#subst#',subst_str); % replace #subst#
    end
    targetReacs = [targetReacs new_target]; % add to target regions
end
targetReacs = targetReacs(setdiff(1:length(targetReacs),target_w_subst_placeholder)); % remove target regions with #subst# placeholder
% replace #subst# in desired regions
desired_w_subst_placeholder = find(cellfun(@(x) any(any(contains(x,'#subst#'))),desiredReacs));
for i = desired_w_subst_placeholder % if one region definition contains #subst#, replace it
    subst_str = ['(' strjoin({subst(:).fac_and_name},' ') ')'];
    desiredReacs{i} = strrep(desiredReacs{i},'#subst#',subst_str); % replace #subst#
end
% OLD AND PROBABLY INCORRECT:
% If desired regions are split up, this would mean that ALL new Desired regions must be
% maintained. But in fact it is sufficient if ONE Desired region remains feasible.
% for i = desired_w_subst_placeholder % if one region definition contains #subst#, replace it
%     comb = cellfun(@(x) double(isempty(x)):1, {subst(:).subst_export_constr},'UniformOutput',0);
%     comb = combvec(comb{:}); % vector that determines which reactions are active in which target region
%     comb = comb(:,any(comb,1)); % remove empty substrate sets
%     new_desired = repmat(desiredReacs(i),1,size(comb,2));
%     for j = 1:size(comb,2)
%         subst_idx = find(comb(:,j));
%         new_desired{j} = [new_desired{j}; vertcat(subst(subst_idx).subst_export_constr)]; % add export constraints
%         subst_str = ['(' strjoin({subst(subst_idx).fac_and_name},' ') ')']; % build string for substrate combination
%         new_desired{j} = strrep(new_desired{j},'#subst#',subst_str); % replace #subst#
%     end
%     desiredReacs = [desiredReacs new_desired]; % add to desired regions
% end
% desiredReacs = desiredReacs(setdiff(1:length(desiredReacs),desired_w_subst_placeholder)); % remove target regions with #subst# placeholder

for i = 1:length(targetReacs)
    [T{i},t{i}] = genV(targetReacs{i},cellstr(cnap.reacID),cnap);
end

% generate d, D
for i = 1:length(desiredReacs)
    [D{i},d{i}] = genV(desiredReacs{i},cellstr(cnap.reacID),cnap);
end

% export idx
target = readOut(table_prod,'targetR',1);
target = strtrim(strsplit(target{1},'/'));
fac_and_prod = strsplit(target{1});
idx.prod = findStrPos(cnap.reacID,fac_and_prod(end));
for i = 1:size(subst,2)
    idx.subs(i) = findStrPos(cnap.reacID,readOut(table_subst{i},'substrate',1));
    idx.subsYieldFactor(i) = str2double(readOut(table_subst{i},'yield_factor',1));
end

% make shure product export is open
if reacMax(idx.prod) == 0
    warning('Product export had upper bound of zero. Was set to 1000.');
    reacMax(idx.prod) = 1000;
end
end

function C = readOut(table,keyword,cols) % returns all cells below the given keyword until the first empty line
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),keyword)));
    for i = 1:length(sheet)
        lastrow = row(i)+find(strcmp(strtrim(table(row(i):end,col(i),sheet)),''),1,'first')-2;
        if isempty(lastrow)
            lastrow = size(table,1);
        end
        C{i} = table((row(i)+1):lastrow,col(i):(col(i)+cols-1),sheet(i));
    end
    if ~exist('C','var')
        error(['ERROR: Keyword ''' keyword ''' was not found in any sheet.'])
    elseif length(C) == 1
        C = C{:};
    else
        disp(['WARNING: Keyword ''' keyword ''' was found multiple times in the sheets.']);
    end
end


function [V,v] = genV(constraints,reacID,cnap)
% Generate Vectors V and v so that V*r <= v
% input: some constraint seperated 3 three cells. e.g.:
%           r_1 + r_4 / r_3 - r_2    |    >=    |    a

    V = [];
    v = [];
    rMin = nan(cnap.numr,1);
    rMax = nan(cnap.numr,1);

    for j = 1:size(constraints,1)
        % get right hand side
        a = constraints{j,3};
        if ischar(a)
            a = str2double(a);
            if isnan(a)
                error('error in ''t'' of target region definition');
            end
        end
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
        numDiv = strtrim(split(constraints{j,1},'/'));
        % get variables and coefficients for vector
        [num, cnum] = findReacAndCoeff(numDiv(1),reacID);
            % fractional constraint ispreprocessed
        if length(numDiv) == 2
            [div, cdiv] = findReacAndCoeff(numDiv(2),reacID);
            [rMin(div), rMax(div)] = CNAfluxVariability(cnap,[],[],-1,div,[],[],0);
            % check if all reactions take identical signs (+ or -)
            % this is needed to do the equation rearrangement. If the signs
            % are ambigous, a lot of case differentiations would be needed,
            % what is not done here.
            if any(rMin(div).*rMax(div) < 0)
                error(['reactions that are part of the divisor inside the '... 
                        'target constraints must not span a positive AND negative range']);
            end
            rDir = sign(rMin(div)+rMax(div));
            if any(rDir.*cdiv' > 0) &&  any(rDir.*cdiv' < 0)
                error(['reactions that are part of the divisor inside the '...
                       'target constraints must all have the same direction '...
                       '(positive or negative). Please check if coefficients and '...
                       'reaction ranges lead to all positive or all negative variables.']);
            end
            if sign(sum(rDir.*cdiv')) == -1 % if the divisor is all negative
                % change direction of inequality
                eqop = -eqop;
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
end

function [ridx,coeff] = findReacAndCoeff(eq,reacID)
    coeff = [];
    r = cell.empty(1,0);
    ridx = [];
    for strPart = strsplit(char(eq))
        str = char(regexprep(strPart,'^(\s|\d|-|\.|\()*|(\s|\d|-|\.|\))*$','')); % remove leading and tailing special characters
        if ~isempty(str)
            v = regexp(reacID, ['^' str '$'], 'match');
            if any(~cellfun(@isempty,v))
                r(end+1) = {str};
                ridx = [ridx find(~cellfun(@isempty,v))'];
            end
        end
    end
    for k = 1:length(r)
        c = regexp(eq, ['(\s|\d|-|\.)*?(?=' r{k} ')'], 'match');
        c = regexprep(char(c{:}),'\s','');
        switch c
            case ''
                coeff(k) = 1;
            case '+'
                coeff(k) = 1;
            case '-'
                coeff(k) = -1;
            otherwise
                coeff(k) = str2double(c);
        end
    end
end