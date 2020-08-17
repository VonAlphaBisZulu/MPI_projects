startup_CNA_parpool;
cnap = CNAloadNetwork({'../network_dirs/iML1515';1},1,1);

options.milp_solver          = 'matlab_cplex'; % 'matlab_cplex'; 
% options.milp_time_limit      = inf;  % 14400; % 4 Stunden ; 72000 = 20 Stunden; 39600 = 11 Stunden; 200000 2.5 Tage
% options.mcs_search_mode      = 1;
% max_solutions                = 5;
% max_num_interv               = 200;
% 
% prod_id = 75;
% subs_id = 1;
% subs_id = [01,14,19];
filepath = './StrainBooster/my_mcs_setups/';

[rmcs, cnap, gmcs, gcnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, file_prod, file_subs] = ...
    computeMCS_fromXLS(cnap,prod_id,subs_id,filepath,max_solutions,max_num_interv,options);