function [rmcs, cnap, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx_cmp_full, prod, subs] = computeMCS_fromXLS(cnap,prod_id,subs_id,filepath,max_solutions,max_num_interv,options);

%% 1. Identifying product and substrate files
prod_id = num2str(prod_id,'%02i');
subs_id = arrayfun(@(x) num2str(x,'%02i'),subs_id,'UniformOutput', false);

filesinPath = dir(filepath);
prod = regexp({filesinPath.name}, ['^P',prod_id,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
subs = regexp({filesinPath.name}, strjoin(strcat('^S',subs_id,'_.*.xls.*'), '|'), 'match'); % SXX_ABC
subs = [subs{~cellfun(@isempty,subs)}];
[~,prod_name] = fileparts(prod{:});
[~,model_name] = fileparts(cnap.path);
filename=['StrainBooster/my_mcs_results/' model_name '-gMCS-' prod_name '-' datestr(date,'yyyy-mm-dd')];

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

cnap.mcs.T = T;
cnap.mcs.t = t;
cnap.mcs.D = D;
cnap.mcs.d = d;
cnap.mcs.rkiCost = rkiCost;
cnap.mcs.rkoCost = rkoCost;
cnap.mcs.gpr_rules = gpr_rules;
cnap.mcs.genes = genes;
cnap.mcs.gkiCost = gkiCost;
cnap.mcs.gkoCost = gkoCost;
cnap.mcs.max_num_interv = max_num_interv;
cnap.mcs.options = options;

%% This can be used to compute the target yield (per biomass or glucose uptake)
idx.bm = find(~~cnap.objFunc);
idx.prodYieldFactor = T{:}(1,idx.prod);
% fv = CNAoptimizeFlux(cnap);
% fv(~cnap.objFunc) = nan;
% fv = 0.3*fv;
% Ymax = CNAoptimizeYield(cnap,full(sparse(1,idx.prod,T{:}(1,idx.prod),1,cnap.numr)),full(sparse(1,idx.bm,1,1,cnap.numr)),fv);
% Y_thresh = Ymax * 0.3;
% disp(Y_thresh);

%% 6. Compute MCS
[rmcs, gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx_cmp_full] = ...
    CNAgeneMCSEnumerator2(cnap,T,t,D,d,rkoCost,rkiCost,max_solutions,max_num_interv,gkoCost,gkiCost,gpr_rules,options,1);

cnap.mcs.rmcs = rmcs;

[~,gselection] = unique(gmcs_idx);
gmcs_selection = gmcs(:,gselection);
gko_text = cell(1,size(gmcs_selection,2));
for i = 1:size(gmcs_selection,2)
    kis = find(~isnan(gmcs_selection(:,i)) & gmcs_selection(:,i) > 0);
    kos = find(~isnan(gmcs_selection(:,i)) & gmcs_selection(:,i) < 0);
    for j = 1:length(kis)
        gko_text{j,i} = ['+ ' strtrim(gcnap.reacID(kis(j),:))];
    end
    gko_text{length(kis)+1,i} = '';
    for j = length(kis)+2:length(kis)+length(kos)
        gko_text{j,i} = ['- ' strtrim(gcnap.reacID(kos(j-length(kis)-1),:))];
    end
end
disp('some gene MCS:');
disp(gko_text);

gmcs = sparse(gmcs);
save([filename '.mat'],'cnap', 'rmcs', 'gmcs', 'gcnap', 'cmp_gmcs', 'cmp_gcnap', 'mcs_idx_cmp_full', 'gko_text', 'gpr_rules','-v7.3');

%% 7. Characterization and ranking of MCS
% Instead of the gene-MCS, their corresponding reaction-representations are analyzed.
% This is preferred, because the reaction-model is smaller and therefore analysis is 
% faster than in the GPR-extended model. Furthermore different gene-MCS can lead to 
% identical 'phenotypes' when translated to the reaction-model and by analyzing rMCS
% only a reduced, non-redundant set of reaction-MCS needs therefore to be considered.
if full(~all(all(isnan(gmcs)))) % if mcs have been found
    disp('Characterizing mcs');
  % 5.1) Lump redundant MCS and create flux bounds for each mutant model
    rmcs(isnan(rmcs)) = -inf; % this step allows to apply 'unique' too remove duplicates
    [rmcs,~,gmcs_rmcs_map] = unique(rmcs','rows');
    rmcs = rmcs';
    rmcs(rmcs == -inf) = nan;
    MCS_mut_lb = repmat({cnap.reacMin},1,size(rmcs,2));
    MCS_mut_ub = repmat({cnap.reacMax},1,size(rmcs,2));
    MCS_mut_lb = arrayfun(@(x) MCS_mut_lb{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);
    MCS_mut_ub = arrayfun(@(x) MCS_mut_ub{x}.*(rmcs(:,x)==1 | rmcs(:,x)==0),1:size(rmcs,2),'UniformOutput',0);
  % 5.2) Set relevant indices [criterion 2-7] and prepare thermodynamic (MDF) parameters [criterion 9]
    % reaction indices
    [idx,mdfParam] = relev_indc_and_mdf_Param(cnap,idx);

  % 5.3) Define core metabolism [criterion 8]
    % Add the new reactions also to the list of reactions that will be
    % considered "core" reactions in the final MCS characterization and ranking
    new_reacs = ismember(cellstr(cnap.reacID),{'ACLDC','BTDD','EX_23bdo_e'});
    reac_in_core_metabolism = cell2mat(CNAgetGenericReactionData_as_array(cnap,'core_reac'));
    reac_in_core_metabolism(isnan(reac_in_core_metabolism)) = 1;
    reac_in_core_metabolism = logical(reac_in_core_metabolism);
    reac_in_core_metabolism(new_reacs) = 1;
    lbCore = cnap.reacMin;
    ubCore = cnap.reacMax;
    lbCore(~reac_in_core_metabolism) = 0;
    ubCore(~reac_in_core_metabolism) = 0;
  % 5.4) Costs for genetic interventions  [criterion 10]
    intvCost                  = gcnap.mcs.kiCost;
    intvCost(isnan(intvCost)) = gcnap.mcs.koCost(isnan(intvCost));
    intvCost(gcnap.rType == 'g') = 1;
    gene_and_reac_names = cellstr(gcnap.reacID);
    gene_and_reac_names(gcnap.rType == 'g') = cellstr(gcnap.reacID((gcnap.rType == 'g'),4:end)); % to avoid the 'GR-' prefix
    gene_and_reac_names(gcnap.rType == 'g') = strrep(strrep(gene_and_reac_names(gcnap.rType == 'g'),'(',''),')',''); % remove brackets
  % 5.5) Start characterization and ranking
    [MCS_rankingStruct, MCS_rankingTable]...
        = CNAcharacterizeGeneMCS( cnap , MCS_mut_lb, MCS_mut_ub, 1:size(MCS_mut_lb,2),... model, mutants LB,UB, incices ranked mcs
        idx, idx.cytMet, D, d, T, t, mdfParam, ... relevant indices, Desired and Target regions
        lbCore, ubCore, gmcs, intvCost, gene_and_reac_names, gmcs_rmcs_map, ...
        0:10, ones(1,10),2); % assessed criteria and weighting factors
    % save ranking and textual gmcs as tab-separated-values
    cell2csv([filename '.tsv'],MCS_rankingTable,char(9));
    text_gmcs = cell(size(gmcs,2),1);
    for i = 1:size(gmcs,2)
        kos = find(~isnan(gmcs(:,i)) & gmcs(:,i) ~= 0);
        for j = 1:length(kos)
            text_gmcs(i,j) = cellstr(gcnap.reacID(kos(j),:));
        end
    end
    cell2csv([filename '-gmcs.tsv'],text_gmcs,char(9));
    save([filename '.mat'],'MCS_rankingStruct','MCS_rankingTable','-append');
end
end

function [idx,mdfParam] = relev_indc_and_mdf_Param(cnap,idx)
% function is used to find reaction and species indices that are used for
% the characterization and ranking of MCS
    % relevant reaction indices
    idx.o2      = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*EX_o2_e.*','match')));
    idx.atpm    = find(~cellfun(@isempty ,regexp(cellstr(cnap.reacID),'.*ATPM','match')));
    % relevant species indices
    idx.cytMet = find(~cellfun(@isempty ,regexp(cellstr(cnap.specID),'_c$','match')));
    % other important species
%     idx.pi    = find(strcmp(cellstr(cnap.specID), 'pi_c')); 
% (adding idx.pi activates the output of coupling mechanism analysis (not yet functional))
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
    % MDF setup (thermodynamic benchmark)
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
    mdfParam.G0 = cell2mat(CNAgetGenericReactionData_as_array(cnap,'deltaGR_0'));
    mdfParam.uncert = cell2mat(CNAgetGenericReactionData_as_array(cnap,'uncertGR_0'));
end
