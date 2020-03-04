e = full( sparse(1,[ findStrPos(cmp_gcnap.reacID,'EX_ac_2_e') ...
                     findStrPos(cmp_gcnap.reacID,'EX_glc__D_2_e') ...
                     findStrPos(cmp_gcnap.reacID,'EX_glyc_2_e') ],[-2,-6,-3],1,cmp_gcnap.numr));
b = full( sparse(1, findStrPos(cmp_gcnap.reacID,'BIOMASS.*','regex'), 1 ,1,cmp_gcnap.numr));
p = full( sparse(1, findStrPos(cmp_gcnap.reacID,'EX_14bdo_e'), 4 ,1,cmp_gcnap.numr));
char(cellfun(@(x) strjoin(cellstr(gcnap.reacID(~isnan(x) & x~=0,:)),char(9)),num2cell(gmcs,1),'UniformOutput',false)');

for i = 1:size(cmp_gmcs,2)
    ko(i) = sum(cmp_gcnap.mcs.koCost(cmp_gmcs(:,i) == -1));
    ki(i) = sum(cmp_gcnap.mcs.kiCost(cmp_gmcs(:,i) == 1));
    cost(i) = ko(i)+ki(i);
end

for i = 1:size(cmp_gmcs,2)
    all_kos = [num2cell(find(isnan(cmp_gmcs(:,i)) | cmp_gmcs(:,i) == -1)) cellstr(cmp_gcnap.reacID(isnan(cmp_gmcs(:,i)) | cmp_gmcs(:,i) == -1,:))];
    kos = find(isnan(cmp_gmcs(:,i)) | cmp_gmcs(:,i) == -1);
    steps = 12;
    index = 1;
    for j = 1:size(kos,2)
        l = kos(:,j);
        l = l(l~=0);
        clear('y');
        cnap_temp = cmp_gcnap;
        cnap_temp.reacMin(l) = 0;
        cnap_temp.reacMax(l) = 0;
        testRegionFeas(cnap_temp,[],cnap_temp.mcs.T,cnap_temp.mcs.t,cnap_temp.mcs.D,cnap_temp.mcs.d);
        ys = CNAplotYieldSpace(cnap_temp,b,e,p,e,steps,[],[],2);
        flipvec = [1:2:(size(ys,1)/2)*2 flip(2:2:(size(ys,1)/2)*2)];
        close(gcf);
        y = figure('Visible','off');
        fill(ys(flipvec,1),ys(flipvec,2),'r');
        
        ylim(y.Children,[0 0.3]);
        xlim(y.Children,[0 0.0045]);
        set(y,'Color','none');
        set(y.Children,'Color','none');
        set(y.CurrentAxes,'XColor','black');
        set(y.CurrentAxes,'YColor','black');
        
        set(y, 'units', 'normalized');
        Tight = get(y.Children, 'TightInset');
        NewPos = [1.04*Tight(1) 1.03*Tight(2) 0.98*(1-Tight(1)-Tight(3)) 0.98*(1-Tight(2)-Tight(4))];
        set(y.Children, 'Position', NewPos);
        
        saveas(y,['figures/' 'BDO_' num2str(i) '.emf']);
        close(y);
        index = index+1;
    end
end

function testRegionFeas(cnap,c_macro,T,t,D,d)
    for i = 1:length(t)
        if ~isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, T{i}, t{i}))
            disp(['At least one target region (T' num2str(i) ') is infeasible in the original model']);
        end
    end
    for i = 1:length(d)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, D{i}, d{i}))
            disp(['At least one desired region (D' num2str(i) ') is infeasible in the original model']);
        end
    end
end