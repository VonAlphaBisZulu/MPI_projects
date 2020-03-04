if ~exist('cnan','var')
    startcna(1)
end
%% MCS computation parameters

cprodidx = 16;
%csubstidx = [01,14,19];
csubstidx = 14;
max_num_interv  = 15;
max_solutions = 200;
time_limit = 14400; % 14400; % 4 Stunden ; 72000 = 20 Stunden; 39600 = 11 Stunden; 200000 2.5 Tage
cnap = CNAloadNetwork({'iML1515';1},1,1);
% cnap = CNAloadNetwork({'iMLcore';1},1,1);
% load('iJO1366/iJO1366GeneNames.mat');
load('iML1515/iML1515GeneNames.mat');
cnap = replaceGeneNames(cnap,ecoliGeneNames);

%% SLURM parpool
tempdir = getenv('SLURM_TEMP');
if ~isempty(tempdir) && isempty(gcp('nocreate'))
%     pdir = prefdir;
%     dir_last_pos = find(prefdir == '/', 1, 'last');
%     jdir = [pdir(1:dir_last_pos) 'local_cluster_jobs' pdir(dir_last_pos:end) '/sjob-' jobid];
    pdir = [tempdir 'prefs'];
    jdir = [tempdir 'jobs'];
    mkdir(pdir); % locate preferences-directory to tmp path
    mkdir(jdir); % locate job/worker-exchange-directory to tmp path
    [~,mem] = unix('sacct -j $SLURM_JOB_ID --format=reqmem -P -n --noconvert'); % set allocated memory
    setenv('SLURM_REQMEM',strtok(mem,'Mn'));                                    % for the job, as env-variable
    try copyfile([prefdir '/matlabprefs.mat'],pdir), catch, end
    prefdir = pdir; % locate preferences-directory to tmp path
    setenv('SLURM_JOB_STORAGE_LOCATION',jdir);
    disp(['slurm job storage location ''' jdir ''' was created and is used for parpool']);
    clstr = parcluster('local');
    clstr.NumWorkers = feature('numcores');
    clstr.JobStorageLocation = jdir;
    parpool(clstr);
end

if exist('cprodidx','var') && isnumeric(cprodidx) && length(cprodidx) == 1
    cprodidx = num2str(cprodidx,'%02i');
else
    error('Specify product in integer variable cprodidx');
end
if exist('csubstidx','var') && isnumeric(csubstidx)
    csubstidx = arrayfun(@(x) num2str(x,'%02i'),csubstidx,'UniformOutput', false);
else
    error('Specify substrates in integer variable or vector cprodidx');
end

filepath = './_StrainBooster/_My_Simulations/';
filesinPath = dir(filepath);
prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*.xls.*'], 'match'); % PXX_ABC
prod = prod{~cellfun(@isempty,prod)};
subs = regexp({filesinPath.name}, strjoin(strcat('^S',csubstidx,'_.*.xls.*'), '|'), 'match'); % SXX_ABC
subs = [subs{~cellfun(@isempty,subs)}];
[~,prod_name] = fileparts(prod{:});

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

[cnap, ~, genes, gr_rules] = CNAgenerateGERassociation( cnap );

%% Load Cut-Set-Calculation parameters from xls
[T, t, D, d,rkoCost,rkiCost,reacMin,reacMax,gkoCost,gkiCost,idx] = CNAgetgMCScalcParamXls( cnap, prod, subs, genes);
cnap.reacMin = reacMin;
cnap.reacMax = reacMax;

%% MCS filename and flux limit
filename=['_StrainBooster/_My_Simulations/Solutions/' cnap.path '-gMCS-' prod_name '-' datestr(date,'yyyy-mm-dd')];
default_flux_limit = 1000;

[gmcs, gcnap, cmp_gmcs, cmp_gcnap, mcs_idx] = CNAgeneMultiImposedMCSFinder(cnap, T , t , D , d ,...
                                                rkoCost,rkiCost, ... koCost, kiCost
                                                max_num_interv,time_limit,max_solutions,...
                                                1,0, ... use_compression,use_bigM,
                                                1,gkoCost,gkiCost,[],1); % enum_method, gkoCost, gkiCost, gr_rules, verbose

gmcs = sparse(gmcs);
save([filename '.mat'],'cnap', 'gr_rules', 'gcnap', 'gmcs', 'gkoCost', 'gkiCost', 'cmp_gcnap', 'cmp_gmcs', 'mcs_idx', 'T', 't', 'D', 'd','-v7.3');

% rmcs = gmcs2rmcs(gmcs,enzymes,numr);

disp('verifying mcs');
if ~isempty(getenv('SLURM_JOB_ID')), parpool(clstr); end
evalgMCS2;

save([filename '.mat'], 'IS_rankingStruct', 'IS_rankingTable','-append');

parfor i = 1:size(gmcs,2)
    cnap_temp = gcnap;
    cnap_temp.reacMin(isnan(gmcs(:,i)) | gmcs(:,i) == -1) = 0;
    cnap_temp.reacMax(isnan(gmcs(:,i)) | gmcs(:,i) == -1) = 0;
    testRegionFeas(cnap_temp,[], gcnap.mcs.T , gcnap.mcs.t , gcnap.mcs.D , gcnap.mcs.d)
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

function cnap = replaceGeneNames(cnap,ecoliGeneNames)
    grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    for i = ecoliGeneNames'
        grRules = strrep(grRules,i(1),i(2));
    end
    cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',grRules);
end