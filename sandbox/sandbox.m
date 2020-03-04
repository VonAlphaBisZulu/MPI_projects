idx.growth = findStrPos(cnap.reacID,'BIOMASS.*core','regex');
idx.subs   = findStrPos(cnap.reacID,'EX_glc__D_2_e');
idx.prod   = findStrPos(cnap.reacID,'EX_etoh_e');
bestMCS = round(linspace(1,size(gmcs,2),8));

yspace = {};
pplane = {};
for i = bestMCS
    cnap = gcnap;
    cnap.reacMin(isnan(gmcs(:,i)) | gmcs(:,i)==-1) = 0;
    cnap.reacMax(isnan(gmcs(:,i)) | gmcs(:,i)==-1) = 0;

   ys = CNAplotYieldSpace(cnap,iv(cnap.numr,idx.growth)' ,-iv(cnap.numr,idx.subs)',...
                           iv(cnap.numr,idx.prod)'   ,-iv(cnap.numr,idx.subs)',10,[],[],0);
   title(['MCS ' num2str(i)]);
   pp = CNAplotPhasePlane(cnap,[idx.growth, idx.prod],[],[],0,10);
   title(['MCS ' num2str(i)]);
%    xlim([0 0.1]);
%    ylim([0 1]);
   yspace = [yspace {ys}];
   pplane = [pplane {pp}];
end
