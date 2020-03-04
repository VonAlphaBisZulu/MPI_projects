% load('_StrainBooster/_My_Simulations/Mechthild Solutions/iML1515-2019-08-06-glc-glyc-ac/iML1515-gMCS-P14_Isobutanol-2019-08-06.mat');
% load('/scratch/CNA_SVN/_StrainBooster/_My_Simulations/mechthildSolutions/iML1515-2019-07-19-glc-glyc-ac/iML1515-gMCS-P16_2-3-Butanediol-2019-07-19.mat');
load('iML1515/thermo.mat');
%% load according subst and prod files
% cprodidx = 8;
% csubstidx = [01,14,19];
% cprodidx = num2str(cprodidx,'%02i');
% csubstidx = arrayfun(@(x) num2str(x,'%02i'),csubstidx,'UniformOutput', false);
% filepath = './_StrainBooster/_My_Simulations/';
% filesinPath = dir(filepath);
% prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*.xls.*'], 'match'); % PXX_ABC
% prod = prod{~cellfun(@isempty,prod)};
% subs = regexp({filesinPath.name}, strjoin(strcat('^S',csubstidx,'_.*.xls.*'), '|'), 'match'); % SXX_ABC
% subs = [subs{~cellfun(@isempty,subs)}];

if ~exist('cnan','var')
    startcna(1)
end
% parpool(3);
cnap_core = CNAloadNetwork({'iMLcore';1},1,1);

% intvCost
intvCost                  = gcnap.mcs.kiCost;
intvCost(isnan(intvCost)) = gcnap.mcs.koCost(isnan(intvCost));
intvCost(gcnap.rType == 'g') = 1;
gcnap.mcs.intvCost = intvCost;

gene_and_reac_names = cellstr(gcnap.reacID);
gene_and_reac_names(gcnap.rType == 'g') = cnap.genes;

% compute reaction mcs from gene mcs, select  and generate bounds
rmcs = gmcs2rmcs(gmcs,gr_rules,gcnap.mcs.rmap,gcnap.rType);
rmcs(isnan(rmcs)) = -inf;
[rmcs,~,gmcs_rmcs_map] = unique(rmcs','rows');
rmcs = rmcs';
rmcs(rmcs == -inf) = nan;
% generate mutants based on reaction knock outs in the metabolic network
IS_lb = repmat({cnap.reacMin},1,size(rmcs,2));
IS_ub = repmat({cnap.reacMax},1,size(rmcs,2));
IS_lb = arrayfun(@(x) IS_lb{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);
IS_ub = arrayfun(@(x) IS_ub{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);

for i = 1:size(rmcs,2) % check if all mutants are valid (Target infeasible, Desired feasible)
    ccnap = cnap;
    ccnap.reacMin = IS_lb{i};
    ccnap.reacMax = IS_ub{i};
    testRegionFeas(ccnap,[],T,t,D,d);
end

% 1
%prod
table = loadSpecReacXLStoStrArray(char(prod));
[row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),'targetR')));
target = strtrim(strsplit(table{row+1,col,sheet},'/'));
fac_and_prod = strsplit(target{1});
idx.prod = findStrPos(cnap.reacID,fac_and_prod(end));
if length(fac_and_prod) == 2
    idx.prodYieldFactor = str2double(fac_and_prod(1));
else
    idx.prodYieldFactor = 1;
end
added_reacs = readOut(table,'reac_id');
%subs
for i = 1:size(subs,2)
    table = loadSpecReacXLStoStrArray(char(subs(i)));
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),'substrate')));
    idx.subs(i) = findStrPos(cnap.reacID,table(row+1,col,sheet));
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),'yield_factor')));
    idx.subsYieldFactor(i) = str2double(table(row+1,col,sheet));
    added_reacs = [added_reacs; readOut(table,'reac_id')];
end
idx.o2      = findStrPos(cellstr(cnap.reacID),'.*EX_o2_.*','regex');
idx.atpm    = findStrPos(cellstr(cnap.reacID),'.*ATPM','regex');
idx.bm      = findStrPos(cellstr(cnap.reacID),'.*BIOMASS_.*_core_.*','regex');
% 2
cytMet = findStrPos(cellstr(cnap.specID),'_c$','regex');
% 3
idx.pi    = find(strcmp(cellstr(cnap.specID), 'pi'));
idx.h     = find(strcmp(cellstr(cnap.specID), 'h_c'));
idx.h2o   = find(strcmp(cellstr(cnap.specID), 'h2o_c'));
idx.atp   = find(strcmp(cellstr(cnap.specID), 'atp_c'));
idx.adp   = find(strcmp(cellstr(cnap.specID), 'adp_c'));
idx.amp   = find(strcmp(cellstr(cnap.specID), 'amp_c'));
idx.nad   = find(strcmp(cellstr(cnap.specID), 'nad_c'));
idx.nadh  = find(strcmp(cellstr(cnap.specID), 'nadh_c'));
idx.nadp  = find(strcmp(cellstr(cnap.specID), 'nadp_c'));
idx.nadph = find(strcmp(cellstr(cnap.specID), 'nadph_c'));
idx.co2_e = find(strcmp(cellstr(cnap.specID), 'co2_c'));
idx.glc_e = find(strcmp(cellstr(cnap.specID), 'glc__D_e'));
mdfParam.Cmin    = 1e-6*ones(cnap.nums,1);
mdfParam.Cmin(idx.glc_e) = 1e-6;
mdfParam.Cmax = 0.02*ones(cnap.nums,1);
mdfParam.Cmax(idx.co2_e) = 1e-4;
mdfParam.Cmax(idx.glc_e) = 0.055557;
mdfParam.fixed_ratios(1,1:3) = [idx.atp   idx.adp   10];
mdfParam.fixed_ratios(2,1:3) = [idx.adp   idx.amp    1];
mdfParam.fixed_ratios(3,1:3) = [idx.nad   idx.nadh  10];
mdfParam.fixed_ratios(4,1:3) = [idx.nadph idx.nadp  10];
mdfParam.RT = 8.31446*300/1000; % Computation of MDF in kJ
mdfParam.bottlenecks = 0; % change to 1 to compute thermodynamic bottlenecks
ord = findStrPos(cellstr(cnap.reacID),cellstr(thermo(:,1)))';
mdfParam.G0 = nan(cnap.numr,1);
mdfParam.G0 (ord(ord~=0)) = cell2mat(thermo(ord~=0,2));
mdfParam.uncert = nan(cnap.numr,1);
mdfParam.uncert(~isnan(mdfParam.G0)) = 0;
% mdfParam.uncert(ord(ord~=0)) = cell2mat(thermo(ord~=0,3));

% 4
lbCore = cnap.reacMin;
ubCore = cnap.reacMax;
reacID = regexprep(cellstr(cnap.reacID),'_gen(_rev$|$)','');
core = [intersect(reacID,cellstr(cnap_core.reacID)); added_reacs];
notCore = findStrPos(cnap.reacID, setdiff(cellstr(cnap.reacID),core));
lbCore(notCore) = 0;
ubCore(notCore) = 0;

% [IS_rankingStruct, IS_rankingTable]...
% = CNAcharacterizeIS_genes( cnap , IS_lb(1:25), IS_ub(1:25), 1:25, idx, cytMet, D, d, T, t, ...
% mdfParam, lbCore, ubCore, gmcs, gcnap.mcs.intvCost, gene_and_reac_names, gmcs_rmcs_map, [2,3,4,10], ones(1,10),2);
% return;

[IS_rankingStruct, IS_rankingTable]...
= CNAcharacterizeIS_genes( cnap , IS_lb, IS_ub, 1:size(IS_lb,2), idx, cytMet, D, d, T, t, ...
mdfParam, lbCore, ubCore, gmcs, gcnap.mcs.intvCost, gene_and_reac_names, gmcs_rmcs_map, 0:10, ones(1,10),2);

function testRegionFeas(cnap,c_macro,T,t,D,d)
    for i = 1:length(t)
        if ~isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, T{i}, t{i}))
            disp(['At least one target region (T' num2str(i) ') is feasible in the original model']);
        end
    end
    for i = 1:length(d)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, D{i}, d{i}))
            disp(['At least one desired region (D' num2str(i) ') is infeasible in the original model']);
        end
    end
end
function C = readOut(table,keyword) % returns all cells below the given keyword until the first empty line
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),keyword)));
    lastrow = row+find(strcmp(strtrim(table(row:end,col,sheet)),''),1,'first')-2;
    if isempty(lastrow)
        lastrow = size(table,1);
    end
    C = table((row+1):lastrow,col,sheet);
    C = reshape(C,numel(C),1);
end