%% Script to compute MCS using another scaling

%ECC2Comp -1e-1


clear all
%% Adapt MATLAB Path
pathname = fileparts(which('Skript_skal.m'));
rmpath(genpath([pathname '/..']));
addpath(pathname);
addpath([pathname '/..']);

if ~exist('cnan','var')
    startcna(1)
end
if license('test','Distrib_Computing_Toolbox') && isempty(getCurrentTask()) && ...
        (~isempty(ver('parallel'))  || ~isempty(ver('distcomp'))) && isempty(gcp('nocreate')) %#ok<DCRENAME>
    parpool();
    wait(parfevalOnAll(@startcna,0,1)); % startcna on all workers
end


%% ECC2Comp

load('ECC2comp.mat');
clear modules;
clear options;

etoh_rID = find(ismember(cnap.reacID,{'R_EthEx'}));

%Set default values
cnap.reacMin(~isnan(cnap.reacDefault)) = cnap.reacDefault(~isnan(cnap.reacDefault));
cnap.reacMax(~isnan(cnap.reacDefault)) = cnap.reacDefault(~isnan(cnap.reacDefault));

%Target for old MCS function
T = full(sparse(1,etoh_rID,-1 ,1,cnap.numr));
t = -0.5;
maxSize = 3;




%% Computation

% Target for new MCS function
modules{1}.sense = 'target';
modules{1}.type  = 'lin_constraints';
modules{1}.V = T;
modules{1}.v = t;

%Computation with old MCS function and CPLEX with indicator constraints
options = struct;
options.milp_solver = 'matlab_cplex';
options.postproc_verify_mcs = 0;
mcs_old_cpx = CNAMCSEnumerator2(cnap,modules{1}.V,modules{1}.v,[],[],[],[],inf,maxSize,options,0);
disp([num2str(size(mcs_old_cpx,2)) ' MCS with old CPLEX.']);

% Computation with new MCS function and cplex with bigM
disp('Computation with CPLEX');
options = struct;
options.postproc_verify_mcs = 0;
options.milp_solver = 'cplex';
options.milp_bigM =5e3;
options.milp_time_limit = inf;
options.mcs_search_mode = 2;
mcs_new_inlinprog = CNAMCSEnumerator3(cnap,modules,[],[],inf,maxSize,options,1);
disp([num2str(size(mcs_new_inlinprog,2)) ' MCS with new cplex.']);


mcs_new_cpx = mcs_new_inlinprog;

disp('finished.');

disp('Verifying MCS for old function: ');
feasible_old = nan(1,size(mcs_old_cpx,2));
parfor i = 1:size(mcs_old_cpx,2)
    kos = find(logical(mcs_old_cpx(:,i)))';
    V = T;
    v = t;
    for ko = kos
        V = [V; sparse([1,2],[ko,ko],[-1,1],2,cnap.numr)];
        v = [v;0;0];
    end
    feasible_old(i) = testRegionFeas(cnap,V,v,-1);
end
if all(~feasible_old)
    disp('All old MCS are valid');
else
    disp([num2str(sum(feasible_old)) ' old MCS are invalid' ]);
end

disp('Verifying MCS for new function: ');
feasible_new = nan(1,size(mcs_old_cpx,2));
parfor i = 1:size(mcs_new_cpx,2)
    kos = find(logical(mcs_new_cpx(:,i)))';
    V = T;
    v = t;
    for ko = kos
        V = [V; sparse([1,2],[ko,ko],[-1,1],2,cnap.numr)];
        v = [v;0;0];
    end
    feasible_new(i) = testRegionFeas(cnap,V,v,-1);
end
if all(~feasible_old)
    disp('All new MCS are valid');
else
    disp([num2str(sum(feasible_old)) ' new MCS are invalid' ]);
end

% Compare both computations
a1 = compare_mcs_sets(mcs_old_cpx,mcs_new_cpx);
a2 = compare_mcs_sets(mcs_new_cpx,mcs_old_cpx);

disp('======================');
disp('Compare old and new CPLEX:');
disp(['Out of ' num2str(size(mcs_old_cpx,2)) ' possible MCS:']);
disp([num2str(sum(a1==3)) ' are identical.']);
disp([num2str(sum(a2==2)) ' were found (as ' num2str(sum(a1==1)) ' supersets with more cuts).']);
disp([num2str(sum(a2==0)) ' MCS were unique to the ''old'' MCS.']);
disp([num2str(sum(a1==0)) ' MCS were unique to the ''new'' MCS.']);


%Computation with intlinprog and big M
disp('Computation with Intlinprog')
options = struct;
options.postproc_verify_mcs = 0;
options.milp_solver = 'cplex';
options.milp_bigM =5e3;
options.milp_time_limit = inf;
options.mcs_search_mode = 2;
mcs_new_inlinprog = CNAMCSEnumerator3(cnap,modules,[],[],inf,maxSize,options,1);
disp([num2str(size(mcs_new_inlinprog,2)) ' MCS with new Intlinprog.']);
mcs_new_cpx = mcs_new_inlinprog;

disp('finished.');

disp('Verifying MCS for new function: ');
feasible_new = nan(1,size(mcs_old_cpx,2));
parfor i = 1:size(mcs_new_cpx,2)
    kos = find(logical(mcs_new_cpx(:,i)))';
    V = T;
    v = t;
    for ko = kos
        V = [V; sparse([1,2],[ko,ko],[-1,1],2,cnap.numr)];
        v = [v;0;0];
    end
    feasible_new(i) = testRegionFeas(cnap,V,v,-1);
end
if all(~feasible_old)
    disp('All new MCS are valid');
else
    disp([num2str(sum(feasible_old)) ' new MCS are invalid' ]);
end

% Compare Intlinprog Big-M with CPLEX-Indicator

a1 = compare_mcs_sets(mcs_old_cpx,mcs_new_cpx);
a2 = compare_mcs_sets(mcs_new_cpx,mcs_old_cpx);

disp('======================');
disp('Compare old CPLEX and new Intlinprog :');
disp(['Out of ' num2str(size(mcs_old_cpx,2)) ' possible MCS:']);
disp([num2str(sum(a1==3)) ' are identical.']);
disp([num2str(sum(a2==2)) ' were found (as ' num2str(sum(a1==1)) ' supersets with more cuts).']);
disp([num2str(sum(a2==0)) ' MCS were unique to the ''old'' MCS.']);
disp([num2str(sum(a1==0)) ' MCS were unique to the ''new'' MCS.']);