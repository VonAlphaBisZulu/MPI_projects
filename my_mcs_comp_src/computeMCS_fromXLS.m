function [rmcs, cnap, gmcs, gcnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, file_prod, file_subs] = computeMCS_fromXLS(cnap,prod_id,subs_id,filepath,max_solutions,max_num_interv,options);

%% 1. Identifying product and substrate files
prod_id = num2str(prod_id,'%02i');
subs_id = arrayfun(@(x) num2str(x,'%02i'),subs_id,'UniformOutput', false);

filesinPath = dir(filepath);
prod = regexp({filesinPath.name}, ['^P',prod_id,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
subs = regexp({filesinPath.name}, strjoin(strcat('^S',subs_id,'_.*.xls.*'), '|'), 'match'); % SXX_ABC
subs = [subs{~cellfun(@isempty,subs)}];
[~,prod_name] = fileparts(prod{:});
filename=['_StrainBooster/_My_Simulations/Solutions/' cnap.path '-gMCS-' prod_name '-' datestr(date,'yyyy-mm-dd')];

%% 2. Adding species and reactions from files
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
        cprintf([0.8 0.6 0.3],[char(prod) ': no reactions were added to model' newline]);
    end
end

%% 3. Checking mass and and charge balances
check_mass_balance(cnap);

%% 4. Get GPR-associations
[cnap, ~, genes, gpr_rules] = CNAgenerateGPRrules( cnap );

%% 5. Load Cut-Set-Calculation parameters from xls
[T, t, D, d,rkoCost,rkiCost,reacMin,reacMax,gkoCost,gkiCost,idx] = CNAgetgMCScalcParamXls( cnap, prod, subs, genes);
cnap.reacMin = reacMin;
cnap.reacMax = reacMax;

idx.growth = find(~~cnap.objFunc);
fv = CNAoptimizeFlux(cnap);
fv(~cnap.objFunc) = nan;
fv = 0.3*fv;
Ymax = CNAoptimizeYield(cnap,full(sparse(1,idx.prod,T{:}(1,idx.prod),1,cnap.numr)),full(sparse(1,idx.growth,1,1,cnap.numr)),fv);
Y_thresh = Ymax * 0.3;
disp(Y_thresh);

%% 6. Compute MCS
[rmcs, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx_cmp_full] = ...
    CNAgeneMCSEnumerator2(cnap,T,t,D,d,rkoCost,rkiCost,max_solutions,max_num_interv,gkoCost,gkiCost,gpr_rules,options,1);

gmcs = sparse(gmcs);
save([filename '.mat'],'cnap', 'rmcs', 'gmcs', 'gcnap', 'cmp_gmcs', 'cmp_gcnap', 'mcs_idx_cmp_full', 'gpr_rules','-v7.3');

evalgMCS2;

save([filename '.mat'], 'IS_rankingStruct', 'IS_rankingTable','-append');

end


