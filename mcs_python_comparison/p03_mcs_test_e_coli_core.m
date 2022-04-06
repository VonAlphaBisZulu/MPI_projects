cnap = CNAsbmlModel2MFNetwork('C:\Users\Philipp\Documents\Python\mcs\examples\ECC2_23BDO.sbml');
load('C:\Users\Philipp\Documents\Python\mcs\matlab\iJO1366geneNames.mat')
cnap.reacID = cnap.reacID(:,3:end);

gene_mcs = false;
special_kos = false;

modules{1}.type = 'lin_constraints';
modules{1}.sense = 'desired';
modules{1}.V = -ismember(cellstr(cnap.reacID),'BIOMASS_Ec_iJO1366_core_53p95M')';
modules{1}.v = -0.05;
modules{2}.type = 'lin_constraints';
modules{2}.sense = 'target';
modules{2}.V = ismember(cellstr(cnap.reacID),'EX_23bdo_e')';
modules{2}.V(ismember(cellstr(cnap.reacID),'EX_glc__D_e')) = 0.3;
modules{2}.v = 0;
maxSolutions = inf;
options.milp_bigM = 0;
options.milp_time_limit = inf;
options.milp_solver = 'gurobi';
options.mcs_search_mode = 3;

if gene_mcs
    for i = 1:length(ecoliGeneNames)
        cnap.reacNotes = strrep(cnap.reacNotes,ecoliGeneNames(i,1),ecoliGeneNames(i,2));
    end
    [~,~,genes] = CNAgenerateGPRrules(cnap);
end

if special_kos
    if gene_mcs
        maxCost = 10;
        koCost = nan(cnap.numr,1);
        koCost(ismember(cellstr(cnap.reacID),'EX_o2_e'	)) = 0.4;
        kiCost = nan(cnap.numr,1);
        kiCost(ismember(cellstr(cnap.reacID),'EX_ac_e'	)) = 0.4;
        gkoCost = nan(numel(genes),1);
        gkoCost(ismember(genes,'aceB'     )) = 1;
        gkiCost = nan(numel(genes),1);
        gkiCost(ismember(genes,'crr'     )) = 1;
    else
        maxCost = 10;
        koCost = nan(cnap.numr,1);
        koCost(ismember(cellstr(cnap.reacID),'PFK'    )) =  1;
        koCost(ismember(cellstr(cnap.reacID),'PFL'    )) =  2;
        koCost(ismember(cellstr(cnap.reacID),'PGI'    )) =  1.4;
        koCost(ismember(cellstr(cnap.reacID),'PGK'    )) =  1.2;
        koCost(ismember(cellstr(cnap.reacID),'PGL'    )) =  2.6;
        koCost(ismember(cellstr(cnap.reacID),'AKGDH'  )) =  1.3;
        koCost(ismember(cellstr(cnap.reacID),'ATPS4r' )) =  0.7;
        koCost(ismember(cellstr(cnap.reacID),'PTAr'   )) =  0.9;
        koCost(ismember(cellstr(cnap.reacID),'PYK'    )) =  0.4;
        koCost(ismember(cellstr(cnap.reacID),'SUCCt3' )) =  1.2;
        koCost(ismember(cellstr(cnap.reacID),'ETOHt2r')) =  2.1;
        koCost(ismember(cellstr(cnap.reacID),'SUCDi'  )) =  0.24;
        koCost(ismember(cellstr(cnap.reacID),'SUCOAS' )) =  0.77;
        koCost(ismember(cellstr(cnap.reacID),'TALA'   )) =  1.5;
        koCost(ismember(cellstr(cnap.reacID),'EX_o2_e')) =  0.9;

        kiCost = nan(cnap.numr,1);
        kiCost(ismember(cellstr(cnap.reacID),'ACALD'   )) = 0.5;
        kiCost(ismember(cellstr(cnap.reacID),'AKGt2r'  )) = 0.4;
        kiCost(ismember(cellstr(cnap.reacID),'PGM'     )) = 0.6;
        kiCost(ismember(cellstr(cnap.reacID),'PIt2r'   )) = 1.3;
        kiCost(ismember(cellstr(cnap.reacID),'ALCD2x'  )) = 0.7;
        kiCost(ismember(cellstr(cnap.reacID),'ACALDt'  )) = 1.8;
        kiCost(ismember(cellstr(cnap.reacID),'ACKr'    )) = 0.12;
        kiCost(ismember(cellstr(cnap.reacID),'PPC'     )) = 0.13;
        kiCost(ismember(cellstr(cnap.reacID),'CS'      )) = 0.2;
        kiCost(ismember(cellstr(cnap.reacID),'RPI'     )) = 0.1;
        kiCost(ismember(cellstr(cnap.reacID),'SUCCt2_2')) = 0.6; 
        kiCost(ismember(cellstr(cnap.reacID),'CYTBD'   )) = 1.3; 
        kiCost(ismember(cellstr(cnap.reacID),'D_LACt2' )) = 0.7; 
        kiCost(ismember(cellstr(cnap.reacID),'ENO'     )) = 1.8;
        kiCost(ismember(cellstr(cnap.reacID),'THD2'    )) = 0.12;
        kiCost(ismember(cellstr(cnap.reacID),'TKT1'    )) = 0.13;
        kiCost(ismember(cellstr(cnap.reacID),'O2t'     )) = 0.2;
        kiCost(ismember(cellstr(cnap.reacID),'PDH'     )) = 0.1;
    end
else
    if gene_mcs
        maxCost = 4;
        koCost = nan(cnap.numr,1);
        koCost(ismember(cellstr(cnap.reacID),'EX_o2_e')) = 0.4;
        kiCost = nan(cnap.numr,1);
        gkoCost = ones(numel(genes),1);
        gkiCost = nan(numel(genes),1);
    else
        maxCost = 4;
        koCost = ones(cnap.numr,1);
        kiCost = nan(cnap.numr,1);
    end
end

if gene_mcs
    [rmcs, full_mcs, full_cnap] = CNAgeneMCSEnumerator3(cnap,modules,...
                    koCost,kiCost,...
                    maxSolutions,maxCost,...
                    gkoCost,gkiCost,[],...
                    options);
    cnap = full_cnap;
    koCost = full_cnap.mcs.koCost(:);
    kiCost = full_cnap.mcs.kiCost(:);
    mcs = full_mcs;
else
    [mcs, status] = CNAMCSEnumerator3(cnap,modules,...
                            koCost,kiCost,...
                            maxSolutions,maxCost,...
                            options);
end

for i = 1:size(mcs,2)
    disp(mcs2text(cnap,mcs(:,i),koCost,kiCost))
end

function ko_ki_text = mcs2text(cnap,mcs,koCost,kiCost)
    kos = find(~isnan(mcs) & ~isnan(koCost) & mcs~=0);
    kis = find(~isnan(mcs) & ~isnan(kiCost) & mcs~=0);
    if ~isempty(kis)
        kis = [strjoin(strcat('+',cellstr(cnap.reacID(kis,:))),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos)
        kos = strjoin(strcat('-',cellstr(cnap.reacID(kos,:))),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end