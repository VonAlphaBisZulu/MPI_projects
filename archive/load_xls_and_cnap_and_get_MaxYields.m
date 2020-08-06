%% Collect xls-product files
path = './_StrainBooster/_My_Simulations/';
filesinPath = dir(path);
prod = regexp({filesinPath.name}, '^P[0-6][0-9]_.*(?=.xls)', 'match'); % PXX_ABC
prod = [prod{~cellfun(@isempty,prod)}];
xlsList = strcat(path,prod,'.xls');

%% start cna and initialize model
if ~exist('cnan','var')
    startcna(1);
end
startcna(1);
cnap = CNAloadNetwork({'iML1515';1},1,1);

maxProd = repmat({'',0,0,0,0},0,1);

idx_rsub        = strcmp(strtrim(cellstr(cnap.reacID)), 'EX_glc__D_e');
idx_nadh16pp    = strcmp(strtrim(cellstr(cnap.reacID)), 'NADH16pp');
idx_o2          = strcmp(strtrim(cellstr(cnap.reacID)), 'EX_o2_e');
idx_atpm        = strcmp(strtrim(cellstr(cnap.reacID)), 'ATPM');

%% iterate through xls-files
cnapBU = cnap;
progBar = Progressbar(0,'Testing alternative products',length(prod));

for idx = 1:length(prod)
progBar = progBar.update(idx);
cnap = cnapBU;
% add reactions and species if necessary
try
    cnap = CNAaddSpecsAndReacsFromFile(cnap,char(xlsList(idx)));
    disp([char(prod(idx)) ': Reactions and Species added to the model']);
catch
    disp([char(prod(idx)) ': Nothing was added to the model']);
end

%% find product sink reaction
TableReads = loadSpecReacXLStoStrArray(char(xlsList(idx)));
[trrow,trcol,cMCSSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'targetR')));
productSynthReac = strsplit(char(TableReads(trrow+1,trcol,cMCSSheet)),'/');

if isempty(productSynthReac(1))
    warning(['synthesis reaction' char(prod(idx)) 'was not found']);
end

idx_rprod       = strcmp(strtrim(cellstr(cnap.reacID)), productSynthReac(1));

if isempty(idx_rprod)
    warning(['synthesis reaction' char(prod(idx)) 'was not found in model']);
end

disp([char(prod(idx)) ': Setting boundaries']);
%% set boundaries
[RMinrow,RMinCcol,~] = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmin'));
[RMaxrow,RMaxcol,~]  = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmax'));

% find range in Table
lastrow = RMinrow+find(strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMinrow+find(~strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'last')-1;
end
Rmin = [cellstr(TableReads(RMinrow+1:lastrow,RMinCcol,cMCSSheet)) num2cell(str2double(TableReads(RMinrow+1:lastrow,RMinCcol+1,cMCSSheet)))];

lastrow = RMaxrow+find(strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMaxrow+find(~strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'last')-1;
end
Rmax = [cellstr(TableReads(RMaxrow+1:lastrow,RMaxcol,cMCSSheet)) num2cell(str2double(TableReads(RMaxrow+1:lastrow,RMaxcol+1,cMCSSheet)))];
% set boundaries
for irmn = 1:length(Rmin)
    if ~isempty(strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))
        cnap.reacMin((strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))) = cell2mat(Rmin(irmn,2));
    else
        warning(['Reaction ' Rmin(irmn,1) ' could not me found. Boundaries were not set.']);
    end
end

for irmx = 1:length(Rmax)
    if ~isempty(strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))
        cnap.reacMax((strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))) = cell2mat(Rmax(irmx,2));
    else
        warning(['Reaction ' Rmax(irmx,1) ' could not me found. Boundaries were not set.']);
    end
end

if cnap.reacMax(idx_rprod) == 0
    warning(['sink reaction for ' char(prod(idx)) ' was 0. Was set to 1000']);
    cnap.reacMax(idx_rprod) = 1000;
end

% prepare FBAs

product             = zeros(1,size(cnap.reacID,1));
product(idx_rprod)  = 1; % positive = yield maximization, negative = minimization
substrate           = zeros(1,size(cnap.reacID,1));
substrate(idx_rsub) = -1;

fixedFluxes                 = zeros(1,size(cnap.reacID,1));
fixedFluxes(idx_rsub)       = -10;
fixedFluxes(fixedFluxes==0) = NaN;
fixedFluxes(idx_atpm) = 0;

maxProd(idx,1) = prod(idx);

% [maxyield,flux_vec,~, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
% disp(['Max yield of ' char(prod(idx)) ' aerobic with reversible NADH16pp: ' num2str(maxyield)])
% disp(['NADH16pp reaction rate: ' num2str(flux_vec(idx_nadh16pp))]);
% maxProd(idx,2) = {maxyield};
% maxProd(idx,6) = {flux_vec(idx_nadh16pp)};

cnap.reacMin(idx_nadh16pp) = 0;
[maxyield,flux_vec,~, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
disp(['Max yield of ' char(prod(idx)) ' aerobic with irreversible NADH16pp: ' num2str(maxyield)]);
%disp(['NADH16pp reaction rate: ' num2str(flux_vec(idx_nadh16pp))]);
%disp(['ATPM rate: ' num2str(flux_vec(idx_atpm))]);
maxProd(idx,3) = {maxyield};
maxProd(idx,7) = {flux_vec(idx_nadh16pp)};
% 
% cnap.reacMin(idx_o2) = 0;
% cnap.reacMin(idx_nadh16pp) = -1000;
% [maxyield,flux_vec,~, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
% disp(['Max yield of ' char(prod(idx)) ' anaerobic with reversible NADH16pp: ' num2str(maxyield)])
% disp(['NADH16pp reaction rate: ' num2str(flux_vec(idx_nadh16pp))]);
% maxProd(idx,4) = {maxyield};
% maxProd(idx,8) = {flux_vec(idx_nadh16pp)};
% cnap.reacMin(idx_nadh16pp) = 0;
% 
% cnap.reacMin(idx_nadh16pp) = 0;
% [maxyield,flux_vec,~, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
% disp(['Max yield of ' char(prod(idx)) ' anaerobic with reversible NADH16pp: ' num2str(maxyield)])
% disp(['NADH16pp reaction rate: ' num2str(flux_vec(idx_nadh16pp))]);
% maxProd(idx,5) = {maxyield};
% maxProd(idx,9) = {flux_vec(idx_nadh16pp)};
% cnap.reacMin(idx_nadh16pp) = 0;
end
nMaxProd = cell2mat(maxProd(:,2:9));
cnap = cnapBU;
progBar.delete();
