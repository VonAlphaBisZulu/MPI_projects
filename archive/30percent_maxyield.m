%% July 2019
% script to identify maximum carbon yield for a specific product on glucose
% return 30% of maximum yield. This threshold is used for later MCS computation
if ~exist('cnapBU','var')
    cnapBU = CNAloadNetwork({'iML1515';1},1,1);
end
parfor i = 1:65
    cprodidx = num2str(i,'%02i');
    c_prod = 1;
    filepath = './_StrainBooster/_My_Simulations/';
    filesinPath = dir(filepath);
    prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*.xls.*'], 'match'); % PXX_ABC
    prod = prod{~cellfun(@isempty,prod)};

    %% network
    cnap = cnapBU;
    try
        cnap = CNAaddSpecsAndReacsFromFile(cnap,prod{:});
    catch
    end

    %% product
    table = loadSpecReacXLStoStrArray(prod{:});
    [row,col,sheet]   = ind2sub(size(table),find(strcmp(strtrim(table),'targetR')));
    target = table(row+1,col,sheet);
    target = strtrim(strsplit(target{:},'/'));
    fac_and_prod = strsplit(target{1});
    product_idx = findStrPos(cnap.reacID,fac_and_prod(2));

    %% boundaries
    ex_reacs = find(sum((sum(abs(cnap.stoichMat))==1).*(cnap.stoichMat==-1)));

    specsWithCarbon = regexp(cellstr(cnap.specNotes), '\[.*C([A-Z]|\d).*]', 'match');
    specsWithCarbon = find(~cellfun(@isempty,specsWithCarbon));
    reacsWCarbon  = cnap.reacID(ex_reacs( ismember(ex_reacs,find(sum(cnap.stoichMat(specsWithCarbon,:))))),:);

    % RMin

    cnap.reacMin(findStrPos(cnap.reacID,reacsWCarbon)) = 0;
    cnap.reacMin(findStrPos(cnap.reacID,'EX_glc__D_e')) = -10;
    cnap.reacMin(findStrPos(cnap.reacID,'EX_co2_e')) = 0; % Under non-uptake conditions
    cnap.reacMin(findStrPos(cnap.reacID,'ATPM')) = 0; % Without ATP-Maintanance

    % RMax
    cnap.reacMax(ex_reacs) = 0;
    cnap.reacMax(findStrPos(cnap.reacID,'BIOMASS.*core','regex')) = 1000;
    cnap.reacMax(findStrPos(cnap.reacID,'BIOMASS.*WT','regex')) = 0;
    cnap.reacMax(findStrPos(cnap.reacID,{   'EX_ac_e'...
                                            'EX_co2_e'...
                                            'EX_etoh_e'...
                                            'EX_for_e'...
                                            'EX_glyc_e'...
                                            'EX_glyc__R_e'...
                                            'EX_h2_e'...
                                            'EX_h2o2_e'...
                                            'EX_h2o_e'...
                                            'EX_h_e'...
                                            'EX_lac__D_e'...
                                            'EX_meoh_e'...
                                            'EX_o2_e'...
                                            'EX_succ_e'...
                                            'EX_tungs_e'})) = 1000;
    cnap.reacMax(findStrPos(cnap.reacID,{   'DM_4crsol_c'...
                                            'DM_5drib_c'...
                                            'DM_aacald_c'...
                                            'DM_amob_c'...
                                            'DM_mththf_c'...
                                            'DM_oxam_c'})) = 0.001;
    c = iv(cnap.numr,product_idx)';
    d = -iv(cnap.numr,findStrPos(cnap.reacID,'EX_glc__D_e'))';
    cnap.reacMax(product_idx) = 1000;

    maxyield(i) = CNAoptimizeYield(cnap,c,d);
    carbyield(i) = maxyield(i) * fac / 6;
    carbyield30percent(i) = 0.3*carbyield(i);
    disp([fac_and_prod{2} ': 30% of max carbon yield is ' num2str(carbyield30percent(i))]);
end