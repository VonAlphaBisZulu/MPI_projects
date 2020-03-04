%load('_StrainBooster/_My_Simulations/ottoSolutions/iJO1366-2StepCalc-P37_L-Lysine-2018-05-31.mat')

idx.growth = findStrPos(cnap.reacID,'BIOMASS_Ec_iJO1366_core_53p95M');
idx.subs   = findStrPos(cnap.reacID,'EX_glc__D_e');
% idx.prod   = findStrPos(cnap.reacID,'EX_14bdo_e');
% bestMCS = [97 240 17 7 146]; % 1,4BDP
 idx.prod   = findStrPos(cnap.reacID,'EX_met__L_e');
 bestMCS = 28; % [42 144 238 244 74]; % Met
%idx.prod   = findStrPos(cnap.reacID,'EX_lys__L_e');
%bestMCS = [35 15 32 33 2]; % Lys
% idx.prod   = findStrPos(cnap.reacID,'EX_glu__L_e');
% bestMCS = [12 43 4 16 26]; % Glu
%mcs = mcs(unique(round(linspace(1,size(mcs,1),1000))),:);

cnapBU = cnap;

yspace = {};
for i = bestMCS
    cnap = cnapBU;
    cnap.reacMin(find(mcs(i,:))) = 0;
    cnap.reacMax(find(mcs(i,:))) = 0;
    cnap.reacMax(idx.subs)  = -1;

   ys = CNAplotYieldSpace(cnap,iv(cnap.numr,idx.growth)' ,-iv(cnap.numr,idx.subs)',...
                           iv(cnap.numr,idx.prod)'   ,-iv(cnap.numr,idx.subs)',10,[],[],2);
   title(['MCS ' num2str(i)])
   xlim([0 0.1]);
   ylim([0 1]);
   yspace = [yspace {ys}];
end
plot_mult

cnap = cnapBU;