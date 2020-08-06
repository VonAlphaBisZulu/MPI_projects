mcss = repmat(~strcmp(cellstr(iJO.reacID),'R_EX_o2_e'),1,size(mcs,1))'.*mcs;
mcso2 = repmat(strcmp(cellstr(iJO.reacID),'R_EX_o2_e'),1,size(mcs,1))'.*mcs;

mcsIndices = round(linspace(1,size(mcs,1),5));

for i = 1:4 %mcsIndices
    disp(iJO.reacID(find(mcss(i,:)),:))
    cprintf('blue',[char(iJO.reacID(find(mcso2(i,:)),:)) char(10) '----' char(10)]);
end