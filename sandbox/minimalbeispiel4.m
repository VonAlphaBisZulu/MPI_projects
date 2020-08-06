cnap = CNAloadNetwork({'iJOcore';1},1,1);

idx.prod = find(strcmp(cellstr(cnap.reacID),'EX_ac_e'));
idx.subs = find(strcmp(cellstr(cnap.reacID),'EX_glc__D_e'));
idx.bm   = find(strcmp(cellstr(cnap.reacID),'BIOMASS_Ec_iJO1366_core_53p95M'));
idx.atpm = find(strcmp(cellstr(cnap.reacID),'ATPM'));
idx.o2   = find(strcmp(cellstr(cnap.reacID),'EX_o2_e'));

T = iv(cnap.numr,idx.prod)'+0.88*iv(cnap.numr,idx.subs)';
t = 0;
D = [-iv(cnap.numr,idx.bm)';-iv(cnap.numr,idx.atpm)'];
d = [-0.05;-3.15];

[cnap, ~] = CNAgenerateGERassociation( cnap );
notknockable = find(~ismember(1:cnap.numr,[cnap.enzymes(:).reactions])); % set reacs w/o gene association to notKOable

maxMCSsize  = 3;
maxMCSnum   = 1e5;
timelimit1 = Inf;
default_flux_limit = 1000;
    
[gmcs,gidx] = CNAgeneMCSEnumerator2(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,[],timelimit1,default_flux_limit,cnap.enzymes);

mcs = CNAgeneToReacMCS(cnap,[],gmcs,gidx);

prod = {'Acetate'};
filepath = './_StrainBooster/_My_Simulations/';

evalgMCS