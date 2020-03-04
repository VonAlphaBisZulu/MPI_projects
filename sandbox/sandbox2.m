% path = '/mechthild/home/schneiderp/Documents/MATLAB/CNA_SVN/_StrainBooster/_My_Simulations/Solutions/iML1515-2019-11-01-GC-glc-glyc-ac/';
% filename = 'iML1515-gMCS-P74_Isobutanol_2NADH_O2-2019-11-01';
% cell2csv([path filename '.xls'],IS_rankingTable,char(9));
ko_text = cell(1,size(gmcs,2));
for i = 1:size(gmcs,2)
    kos = find(~isnan(gmcs(:,i)) & gmcs(:,i) ~= 0);
    for j = 1:length(kos)
        ko_text{j,i} = cellstr(gcnap.reacID(kos(j),:));
    end
end
cmp_ko_text = cell(1,size(cmp_gmcs,2));
for i = 1:size(cmp_gmcs,2)
    kos = find(~isnan(cmp_gmcs(:,i)) & cmp_gmcs(:,i) ~= 0);
    for j = 1:length(kos)
        cmp_ko_text{j,i} = cellstr(cmp_gcnap.reacID(kos(j),:));
    end
end
rmcs = gmcs2rmcs(gmcs,gpr_rules,gcnap.mcs.rmap,gcnap.rType);
rko_text = cell(1,size(rmcs,2));
for i = 1:size(rmcs,2)
    kos = find(~isnan(rmcs(:,i)) & rmcs(:,i) ~= 0);
    for j = 1:length(kos)
        rko_text{j,i} = cellstr(cnap.reacID(kos(j),:));
    end
end