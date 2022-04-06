%% Initialize and load models
if ~exist('cnan','var')
    startcna(1)
end
%% Example 1
% cnap = CNAsbmlModel2MFNetwork('C:\Users\phili\Dokumente\Python\mcs\examples\SmallExample.sbml');
% cnap.reacID = cnap.reacID(:,3:end);
% clear modules
% modules{1}.type = 'lin_constraints';
% modules{1}.sense = 'target';
% modules{1}.V = -ismember(cellstr(cnap.reacID),'R04')';
% modules{1}.v = -1;
% modules{1}.V(2,:) = 2*-ismember(cellstr(cnap.reacID),'R01')';
% modules{1}.v(2,1) = -0.5;
% modules{2}.type = 'lin_constraints';
% modules{2}.sense = 'desired';
% modules{2}.V = -ismember(cellstr(cnap.reacID),'R03')';
% modules{2}.v = -1;
% 
% koCost = [nan nan nan nan 2 3 4 4 4 4]';
% kiCost = [nan 6 nan nan nan nan nan nan nan nan]';
% maxCost = 50;
% maxSolutions = 50;
% options.milp_bigM = 0;
% options.milp_time_limit = inf;
% 
% obj = MCS_cplex(cnap,modules,koCost,kiCost,maxCost,options.milp_bigM);
% [mcs,status] = obj.findMCS(maxSolutions,options.milp_time_limit);

%% Example 2
cnap = CNAsbmlModel2MFNetwork('C:\Users\Philipp\Documents\Python\mcs\examples\iML1515core.sbml');
cnap.reacID = cnap.reacID(:,3:end);
clear modules
modules{1}.type = 'bilev_w_constr';
modules{1}.sense = 'desired';
modules{1}.c = zeros(1,cnap.numr);
modules{1}.c(ismember(cellstr(cnap.reacID),'BIOMASS_Ec_iML1515_core_75p37M')) = -1;
modules{1}.V = -ismember(cellstr(cnap.reacID),'EX_etoh_e')';
modules{1}.v = -1;
modules{1}.V(2,:) = -2*ismember(cellstr(cnap.reacID),'BIOMASS_Ec_iML1515_core_75p37M')';
modules{1}.v(2,1) = -0.1;

koCost = nan(cnap.numr,1);
% koCost = nan(cnap.numr,1);
% koCost(ismember(cellstr(cnap.reacID),{'CYTBO3_4pp','O2Up','h2oEx','AcEx','ATPS4rpp','MDH','PFK','PPS','AcEx','FBA'})) = [2,3,1.5,1.2,1,1,1.1,1,1,1];
kiCost = nan(cnap.numr,1);
% kiCost(ismember(cellstr(cnap.reacID),{'ACKr', 'ICL','PTS','PFL', 'EDD', 'ENO'})) = [2.1,0.9,1.5,3,2,0.7];

koCost(ismember(cellstr(cnap.reacID),{'PTAr' 	})) = 1   ;
koCost(ismember(cellstr(cnap.reacID),{'ACKr'	})) = 2   ;
koCost(ismember(cellstr(cnap.reacID),{'FORtppi'	})) = 1.4 ;
koCost(ismember(cellstr(cnap.reacID),{'SUCDi'	})) = 1.2 ;
koCost(ismember(cellstr(cnap.reacID),{'H2Otpp'	})) = 2.6 ;
koCost(ismember(cellstr(cnap.reacID),{'ATPS4rpp'})) = 1.3 ;
koCost(ismember(cellstr(cnap.reacID),{'PFL'		})) = 0.7 ;
koCost(ismember(cellstr(cnap.reacID),{'GLYCtex'	})) = 0.9 ;
koCost(ismember(cellstr(cnap.reacID),{'EX_ac_e'	})) = 0.4 ;
koCost(ismember(cellstr(cnap.reacID),{'ACtex'   })) = 1.2 ;
koCost(ismember(cellstr(cnap.reacID),{'CYTBO3_4pp'}))=2.1 ;
koCost(ismember(cellstr(cnap.reacID),{'O2tpp'	})) = 0.24;
koCost(ismember(cellstr(cnap.reacID),{'O2tex'	})) = 0.77;
koCost(ismember(cellstr(cnap.reacID),{'EX_o2_e'	})) = 0.4 ;
koCost(ismember(cellstr(cnap.reacID),{'POR5'	})) = 1.5 ;
koCost(ismember(cellstr(cnap.reacID),{'CO2tpp'	})) = 3   ;
koCost(ismember(cellstr(cnap.reacID),{'ALCD19'	})) = 2   ;
koCost(ismember(cellstr(cnap.reacID),{'ASPtpp' 	})) = 0.25;

kiCost(ismember(cellstr(cnap.reacID),{'EDD' 	})) = 0.5 ;
kiCost(ismember(cellstr(cnap.reacID),{'GAPP' 	})) = 0.4 ;
kiCost(ismember(cellstr(cnap.reacID),{'ENO'	 	})) = 0.6 ;
kiCost(ismember(cellstr(cnap.reacID),{'EDA' 	})) = 1.3 ;
kiCost(ismember(cellstr(cnap.reacID),{'AKGDH'  	})) = 0.7 ;
kiCost(ismember(cellstr(cnap.reacID),{'MALD' 	})) = 1.8 ;
kiCost(ismember(cellstr(cnap.reacID),{'NADTRHD' })) = 0.12;
kiCost(ismember(cellstr(cnap.reacID),{'LCARS'	})) = 0.13;
kiCost(ismember(cellstr(cnap.reacID),{'SUCCt1pp'})) = 0.2 ;
kiCost(ismember(cellstr(cnap.reacID),{'SUCFUMtpp'}))=0.1 ;

maxCost = 7;
maxSolutions = inf;
options.milp_bigM = 0;
options.milp_time_limit = inf;
options.solver = 'gurobi';
options.mcs_search_mode = 3;

[mcs_i, status] = CNAMCSEnumerator3(cnap,modules,...
                    koCost,kiCost,...
                    maxSolutions,maxCost,...
                    options);

% obj = MCS_cplex(cnap,modules,koCost,kiCost,maxCost,options.milp_bigM);
% [mcs,status] = obj.findMCS(maxSolutions,options.milp_time_limit);