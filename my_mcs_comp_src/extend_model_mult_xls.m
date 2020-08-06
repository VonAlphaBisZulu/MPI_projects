function [cnap , idx] = extend_model_mult_xls(cnap,cprodidx,csubstidx)
if exist('cprodidx','var') && isnumeric(cprodidx) && length(cprodidx) == 1
    cprodidx = num2str(cprodidx,'%02i');
else
    error('Specify product in integer variable cprodidx');
end
if exist('csubstidx','var') && isnumeric(csubstidx)
    csubstidx = arrayfun(@(x) num2str(x,'%02i'),csubstidx,'UniformOutput', false);
else
    error('Specify substrates in integer variable or vector cprodidx');
end

filepath = './_StrainBooster/_My_Simulations/';
filesinPath = dir(filepath);
prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
subs = regexp({filesinPath.name}, strjoin(strcat('^S',csubstidx,'_.*.xls.*'), '|'), 'match'); % SXX_ABC
subs = [subs{~cellfun(@isempty,subs)}];

disp(['Loading reactions and species from file: ' strjoin(prod,', ')]);
try %% Add new reactions to model from product-xls (if any were defined)
    cnap = CNAaddSpecsAndReacsFromFile(cnap,prod{:});
catch
    cprintf([0.8 0.6 0.3],[char(prod) ': no reactions were added to model' newline]);
end

for file = subs
    disp(['Loading additional substrate uptake pseudoreactions from file: ' (file{:})]);
    try %% Add new reactions to model from substrate-xls (if any were defined)
        cnap = CNAaddSpecsAndReacsFromFile(cnap,file{:});
    catch
        cprintf([0.8 0.6 0.3],[char(file) ': no reactions were added to model' newline]);
    end
end
[~, ~, ~, ~,~,~,cnap.reacMin,cnap.reacMax] = CNAgetgMCScalcParamXls( cnap, prod, subs);
% make sure product export is open
table = loadSpecReacXLStoStrArray(char(prod));
[row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),'targetR')));
target = table(row(1)+1,col(1),sheet(1));
target = strtrim(strsplit(target{:},'/'));
fac_and_prod = strsplit(target{1});
prodct.idx = findStrPos(cnap.reacID,fac_and_prod(2));
prodct.fac = str2double(fac_and_prod(1));

if cnap.reacMax(prodct.idx) == 0
    warning('Product export had upper bound of zero. Was set to 1000.');
    cnap.reacMax(prodct.idx) = 1000;
end

for i = 1:length(subs)
    table_subst = loadSpecReacXLStoStrArray(subs{i});
    substt(i).idx = findStrPos(cnap.reacID,readOut(table_subst,'substrate',1));
    substt(i).fac = str2double(cell2mat(readOut(table_subst,'yield_factor',1)));
end

idx.prod = prodct.idx;
idx.subs = [substt(:).idx];
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

