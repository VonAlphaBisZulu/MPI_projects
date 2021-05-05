function compute_mcs_ext(jdir, filename)
try
    if ~isempty(getenv('SLURM_JOB_ID')) && isempty(gcp('nocreate'))
        tempdir = getenv('SLURM_TEMP');
        if ~isempty(tempdir) && isempty(gcp('nocreate'))
            [~,mem] = unix('sacct -j $SLURM_JOB_ID --format=reqmem -P -n --noconvert'); % set allocated memory
            setenv('SLURM_REQMEM',strtok(mem,'Mn'));                                    % for the job, as env-variable
            try copyfile([prefdir() '/matlabprefs.mat'],pdir), catch, end
            clstr = parcluster('local');
            clstr.NumWorkers = feature('numcores');
            clstr.JobStorageLocation = jdir;
            parpool(clstr);
            wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
        end
    elseif license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
           (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) %#ok<DCRENAME>
        clstr = parcluster('local');
        clstr.JobStorageLocation = jdir;
        parpool();
        wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
    end
    load(filename,'gene_mcs','cnap','modules','rkoCost','maxSolutions','maxCost','gkoCost','options','verbose');
    tic;
    if gene_mcs  %#ok<*UNRCH>
        [rmcs, gmcs, ~, ~, ~, ~, status] = ...
            CNAgeneMCSEnumerator3(cnap,modules,...
                rkoCost,[],...
                maxSolutions,maxCost,...
                gkoCost,[],[],...
                options,verbose);
    else
         [rmcs, status] = CNAMCSEnumerator3(cnap,modules,...
                rkoCost,[],...
                maxSolutions,maxCost,...
                options,verbose);
          gmcs = [];
    end
    comptime = toc;
    save(filename,'gmcs','rmcs','comptime','status');
catch err
    warning('off','backtrace')
    warning( getReport( err, 'extended', 'hyperlinks', 'on' ) ); % getReport not yet implemented in octave
    warning('Something went wrong.');
    warning('on','backtrace')
end
