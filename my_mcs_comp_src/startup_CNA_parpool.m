if ~exist('cnan','var') || isempty(cnan)
    startcna(1);
end
% If runnning on a system with a SLURM workload manager:
% Use directory on internal memory to share data between the workers. 
% If job is running as a SLURM ARRAY, the compression switches (and also other
% parameters if indicated) are overwritten
if ~isempty(getenv('SLURM_ARRAY_TASK_ID')) % overwrite options if a SLURM array is used
    [options,model] = derive_options_from_SLURM_array(str2double(getenv('SLURM_ARRAY_TASK_ID')));
end
if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate'))
    % start parpool and locate preferences-directory to tmp path
    currdir = pwd;
    cd(cnan.cnapath);
    prefdir = start_parallel_pool_on_SLURM_node();
    cd(currdir);
% If running on local machine, start parallel pool and keep compression
% flags as defined above.
%
% On a local machine without a SLURM workload manager, but MATLAB parallel toolbox installed:
elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
       (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate'))
    currdir = pwd;
    cd(cnan.cnapath);
    parpool();
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
    cd(currdir);
end

function pdir = start_parallel_pool_on_SLURM_node()
tempdir = getenv('SLURM_TEMP');
if ~isempty(tempdir) && isempty(gcp('nocreate'))
    pdir = [tempdir 'prefs'];
    jdir = [tempdir 'jobs'];
    mkdir(pdir); % locate preferences-directory to tmp path
    mkdir(jdir); % locate job/worker-exchange-directory to tmp path
    [~,mem] = unix('sacct -j $SLURM_JOB_ID --format=reqmem -P -n --noconvert'); % set allocated memory
    setenv('SLURM_REQMEM',strtok(mem,'Mn'));                                    % for the job, as env-variable
    try copyfile([prefdir() '/matlabprefs.mat'],pdir), catch, end
    setenv('SLURM_JOB_STORAGE_LOCATION',jdir);
    disp(['slurm job storage location ''' jdir ''' was created and is used for parpool']);
    clstr = parcluster('local');
    clstr.NumWorkers = feature('numcores');
    clstr.JobStorageLocation = jdir;
    parpool(clstr);
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end
end