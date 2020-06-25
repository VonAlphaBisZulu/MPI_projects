%% use MATLAB cell mode to successively execute the separate sections
%load XXX.mat

flux2=ECGS;

o2in= find(strcmp(cellstr(flux2.reacID), 'R_EX_o2_LPAREN_e_RPAREN_'));
atpm= find(strcmp(cellstr(flux2.reacID), 'R_ATPM')); % basically ATP -> ADP + P
glucose= find(strcmp(cellstr(flux2.reacID), 'R_EX_glc_LPAREN_e_RPAREN_'));
mue= find(strcmp(cellstr(flux2.reacID), 'R_Ec_biomass_iJO1366_core_53p95M'));
etoh= find(strcmp(cellstr(flux2.reacID), 'R_EX_etoh_LPAREN_e_RPAREN_'));
lactate= find(strcmp(cellstr(flux2.reacID), 'R_EX_lac_DASH_D_LPAREN_e_RPAREN_'));

glucose_uptake_limit= 18.5;
maintenance_limit= 3.15;
product= etoh;
min_yield= 1.4;
min_mue= 0.05;
%no_cuts= [];

n= size(flux2.stoichMat, 2);

inhf= zeros(0, n);
inhf(1, glucose)= -1;
inhf(2, atpm)= -1;
inhf(3, [product, glucose])= [1, min_yield]; % +min_yield because glucose influx has a negative sign!
ub= [glucose_uptake_limit; -maintenance_limit; 0];

%if oxygen_uptake_limit > 0
%  inhf(end+1, o2in)= -1;
%  ub(end+1)= oxygen_uptake_limit;
%end

des= zeros(0, n);
des(1, mue)= -1;
des(2, atpm)= -1;
des(3, glucose)= -1;
des(4, [product, glucose])= -[1, min_yield];
db= [-min_mue; -maintenance_limit; glucose_uptake_limit; 0];

%des=[]
%db=[];

%set finite bounds for useIntegrated==1
zw=find(flux2.reacMin==-inf);
flux2.reacMin(zw)=-10000;
zw=find(flux2.reacMax==inf);
flux2.reacMax(zw)=10000;

% % set all outfluxes to inactive
 flux2.reacMax(1:332)=0;
 % reset reactions from small-scale model to active 
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_Ec_biomass_iJO1366_core_53p95M')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_o2_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_etoh_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_lac_DASH_D_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_ac_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_for_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_succ_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_h2o_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_h2_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_h_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_co2_LPAREN_e_RPAREN_')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_DM_4CRSOL')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_DM_5DRIB')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_DM_AMOB')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_DM_MTHTHF')))=10000;
 flux2.reacMax(find(strcmp(cellstr(flux2.reacID), 'R_EX_meoh_LPAREN_e_RPAREN_')))=10000;

flux2.reacMax(find(flux2.reacMin(1:332)<0))= 10000;
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_EX_succ_LPAREN_e_RPAREN_')))=0;
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_EX_ac_LPAREN_e_RPAREN_')))=  0;
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_EX_glyc_LPAREN_e_RPAREN_')))=0;

% glucose limit to -20 and atp limit to 3.15:
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_EX_glc_LPAREN_e_RPAREN_')))=-20;
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_ATPM')))=3.15;
% R_NADH16pp as reversible
flux2.reacMin(find(strcmp(cellstr(flux2.reacID), 'R_NADH16pp')))= -10000;



flux2.reacMin(o2in)=0;
flux2.reacMax(o2in)=0;
no_cutsGS=[];
           
maxMCSsize=6;
maxMCSnum=1000000;
useIntegrated=0;
filename='EthanolGSMCS040216.txt';

ttt=cputime;
[mcs,reacNames]=CNAregMCSEnumerator(flux2,inhf,ub,des,db,no_cutsGS,maxMCSnum,maxMCSsize,filename,useIntegrated);
%[mcs]=CNAMCSEnumerator(flux2,inhf,ub,des,db,no_cuts,maxMCSnum,maxMCSsize,filename);
cputime-ttt