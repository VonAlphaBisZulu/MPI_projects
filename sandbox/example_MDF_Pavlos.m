%load('/scratch/CNA_SVN/_StrainBooster/_My_Simulations/Solutions/iJO1366-2StepCalc-P15_1-4-Butanediol-2018-06-19.mat')
load('/scratch/CNA_SVN/_StrainBooster/_My_Simulations/Solutions/iJO1366-2StepCalc-P74_L-Methionine-2018-06-19.mat')
idx.atp   = min(findStrPos(cnap.specID,'((?<![^\_])|^)atp((?![^\_])|$)','regex')); % get only first match
idx.adp   = min(findStrPos(cnap.specID,'((?<![^\_])|^)adp((?![^\_])|$)','regex'));
idx.amp   = min(findStrPos(cnap.specID,'((?<![^\_])|^)amp((?![^\_])|$)','regex'));
idx.pi    = min(findStrPos(cnap.specID,'((?<![^\_])|^)pi((?![^\_])|$)','regex'));
idx.h2o   = min(findStrPos(cnap.specID,'((?<![^\_])|^)h2o((?![^\_])|$)','regex'));
idx.h     = min(findStrPos(cnap.specID,'((?<![^\_])|^)h((?![^\_])|$)','regex'));
idx.nad   = min(findStrPos(cnap.specID,'((?<![^\_])|^)nad((?![^\_])|$)','regex'));
idx.nadh  = min(findStrPos(cnap.specID,'((?<![^\_])|^)nadh((?![^\_])|$)','regex'));
idx.nadp  = min(findStrPos(cnap.specID,'((?<![^\_])|^)nadp((?![^\_])|$)','regex'));
idx.nadph = min(findStrPos(cnap.specID,'((?<![^\_])|^)nadph((?![^\_])|$)','regex'));

G0 = cell2mat(CNAgetGenericReactionData_as_array(cnap,'deltaGR_0'));
uncert = cell2mat(CNAgetGenericReactionData_as_array(cnap,'uncertGR_0'));
RT = 8.31446*300/1000; % Computation in kJ
Cmin = 1e-6*ones(cnap.nums,1);
Cmax = 0.02*ones(cnap.nums,1);
Cmax(findStrPos(cnap.specID,'co2_c')) = 1e-4;
Cmin(findStrPos(cnap.specID,'glc__D_e')) = 1e-6;
Cmax(findStrPos(cnap.specID,'glc__D_e')) = 0.055557;
fixed_ratios(1,1:3) = [idx.atp   idx.adp   10];
fixed_ratios(2,1:3) = [idx.adp   idx.amp    1];
fixed_ratios(3,1:3) = [idx.nad   idx.nadh  10];
fixed_ratios(4,1:3) = [idx.nadph idx.nadp  10];
mcs = logical(mcs(10,:));
cnap.reacMin(mcs)=0;
cnap.reacMax(mcs)=0;

[mdf, v, conc, min_df, max_df, dfs,~,~,reac_map]= max_min_driving_force_pathway(cnap.stoichMat,...
  cnap.stoichMat(cnap.specInternal, :), cnap.reacMin, cnap.reacMax, [], [],...
  max(1e-9, cnap.epsilon), [], 1000, RT, [G0 uncert], Cmin, Cmax, fixed_ratios,...
  isnan(G0), false, cnap.reacMin, cnap.reacMax);

[min_df(616) max_df(616) min_df(reac_map==-616) max_df(reac_map==-616)] % warum unterschiedliche Zahlen? Uncertainty?
[min_df(887) max_df(887) min_df(reac_map==-887) max_df(reac_map==-887)]

% min_or_max_r(isnan(dfs))  = {'nan'};
% min_or_max_r(~isnan(dfs)) = {'int!'};
% 
% min_or_max_r(abs(1-abs(dfs./min_df))<0.01) = {'min'};
% min_or_max_r(abs(1-abs(dfs./max_df))<0.01) = {'max'};
% reacs = [cellstr(cnap.reacID(logical(v),:)) num2cell(v(logical(v))) num2cell(dfs(logical(v)))];
% 
% min_or_max_s(logical(exp(conc)')) = {'int!'};
% min_or_max_s(abs(1-abs(exp(conc)'./Cmin))<0.01) = {'min'};
% min_or_max_s(abs(1-abs(exp(conc)'./Cmax))<0.01) = {'max'};
% specs = [cellstr(cnap.specID) num2cell(exp(conc)') min_or_max_s'];