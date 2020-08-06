startcna(1)
prodname = 'etoh';
rid         = ['R_EX_', prodname, '_LPAREN_e_RPAREN_'];

%cnap = stored(strcmp({stored.name},prodname)).iJO;
%mcs  = stored(strcmp({stored.name},prodname)).mcs;


c           = char(cnap.reacID);
clear prody

%% get reaction equations 
    
        % Find educts and products
    %
    %       {'name1'    'educt 1'    'product1'}
    % str = {'name2'    'educt 2'    'product2'}
    %       {'name3'    'educt 3'    'product3'}
    %       {'name4'    'educt 4'    'product4'}
str=[];
for reacIndex = 1:cnap.numr
    educts='';
    zw=find(cnap.stoichMat(:,reacIndex)<0);

    if~isempty(zw)
        educts=[num2str(-cnap.stoichMat(zw(1),reacIndex)),' ',deblank(cnap.specID(zw(1),:))];
        for j=2:length(zw)
            educts=[educts,' + ',num2str(-cnap.stoichMat(zw(j),reacIndex)),' ',deblank(cnap.specID(zw(j),:))];
        end
    end


    products='';
    zw=find(cnap.stoichMat(:,reacIndex)>0);
    if~isempty(zw)
        products=[products,num2str(cnap.stoichMat(zw(1),reacIndex)),' ',deblank(cnap.specID(zw(1),:))];
        for j=2:length(zw)
            products=[products,' + ',num2str(cnap.stoichMat(zw(j),reacIndex)),' ',deblank(cnap.specID(zw(j),:))];
        end
    end

    str=[str; cellstr([educts ' => ' products])];
end

cnap.reacMin(cnap.reacMin<-999) = -Inf;
cnap.reacMax(cnap.reacMax>999) = Inf;
% Save model

cnapbu = cnap;
    
for i = round(linspace(1,size(mcs,1),15))% 1:20 %size(mcs,1)
    disp(i);
    % reload model without cuts
    cnap = cnapbu;
    %% definitions

    subsr = 'R_EX_glc_LPAREN_e_RPAREN_';
    muer  = 'R_Ec_biomass_iJO1366_core_53p95M';
    atpr  = 'R_ATPM';


    %% find reactions and species in model
    idx_rprod       = strcmp(strtrim(cellstr(cnap.reacID)), rid);
    idx_rsub        = strcmp(strtrim(cellstr(cnap.reacID)), subsr);
    idx_ratp        = strcmp(strtrim(cellstr(cnap.reacID)), atpr);
    idx_rmue        = strcmp(strtrim(cellstr(cnap.reacID)), muer);

    %% get knockouts

    imcs             = find(mcs(i,:));
    prody(i).mcs     = imcs;
    

    %% set knocked out reactions to 0
    cnap.reacMin(imcs) = 0;
    cnap.reacMax(imcs) = 0;

    % allow product exchange
    cnap.reacMax(idx_rprod) = Inf;

    % control atp-maintanance
    cnap.reacMin(idx_ratp) = 3.15;
    cnap.reacMax(idx_ratp) = Inf;
    
    % allow substrate uptake
    cnap.reacMin(idx_rsub) = -10;
    cnap.reacMax(idx_rsub) = Inf;

    % biomass rate can optionally be fixed

    cnap.reacMin(idx_rmue) = 0;
    cnap.reacMax(idx_rmue) = Inf;

    %% Define yield optimization
    % yield = product*r/substrate*r

    product             = zeros(1,size(cnap.reacID,1));
    product(idx_rprod)  = 1; % positive = yield maximization, negative = minimization
    substrate           = zeros(1,size(cnap.reacID,1));
    substrate(idx_rsub) = -1;
    muev                = zeros(1,size(cnap.reacID,1));
    muev(idx_rmue)      = 1; % positive = yield maximization, negative = minimization
    atpv                = zeros(1,size(cnap.reacID,1));
    atpv(idx_ratp)      = 1; % positive = yield maximization, negative = minimization
    %product = muev;

    fixedFluxes                 = zeros(1,size(cnap.reacID,1));
    %fixedFluxes(idx_rsub)       = -10;
    fixedFluxes(fixedFluxes==0) = NaN;

    % Y max P/S
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
    prody(i).sup(1) = flux_vec(idx_rsub);
    prody(i).YmaxPS = maxyield;
    
    %% Display results
    output       = [cellstr(cnap.reacID) num2cell(flux_vec) str];
    output_short = output(cell2mat(output(:,2))~=0 & ~isnan(cell2mat(output(:,2))),:);
    output_very_short = output(abs(cell2mat(output(:,2)))>0.05  & ~isnan(cell2mat(output(:,2))),:);

    disp(output_very_short);
    minmax = 'min';
    if sum(product) > 0
        minmax = 'max';
    end
    disp([minmax ' yield = ' num2str(maxyield) ...
        ' (' char(strtrim(cnap.reacID(find(product),:))) ' / ' strtrim(char(cnap.reacID(find(substrate),:))) ')' char(10)...
        num2str(size(output_short,1)) ' active reactions.' char(10)...
        'mue = ' num2str(output{idx_rmue,2}) char(10)...
        'product yield = ' num2str((2*(sum(product)<0)-1)*cell2mat(output(idx_rprod,2))/cell2mat(output(idx_rsub,2)))...
        ' (' strtrim(char(cnap.reacID(find(idx_rprod),:))) ' / ' strtrim(char(cnap.reacID(find(substrate),:))) ')' char(10)...
        strtrim(char(cnap.reacID(find(idx_rprod),:))) ' = '  num2str(cell2mat(output(idx_rprod,2)))]);
    % Y max mue/S
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, muev, substrate,fixedFluxes,[]);
    prody(i).YmaxMueS = maxyield;
    prody(i).sup(2) = flux_vec(idx_rsub);
    disp(['mue max = ' num2str(flux_vec(find(muev)))]);
    
    % Y min P/S
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, -product, substrate,fixedFluxes,[]);
    prody(i).YminPS = -maxyield;
    prody(i).sup(3) = flux_vec(idx_rsub);

    % Y max ATP/S
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, atpv, substrate,fixedFluxes,[]);
    prody(i).YmaxAS = maxyield;
    prody(i).sup(4) = flux_vec(idx_rsub);

    % get min and max ATP-Maintanence rates
    cnap.objFunc(idx_ratp)=1;
    fixedFluxes(idx_rsub) = -10;
    [flux_vec, success, status, optval] = CNAoptimizeFlux(cnap, fixedFluxes);
    prody(i).sup(5) = flux_vec(idx_rsub);
    ratpmin = flux_vec(idx_ratp);

    cnap.objFunc(idx_ratp)=-1;
    [flux_vec, success, status, optval] = CNAoptimizeFlux(cnap, fixedFluxes);
    prody(i).sup(6) = flux_vec(idx_rsub);
    ratpmax = flux_vec(idx_ratp);
    
    prody(i).ratpmin = ratpmin;
    prody(i).ratpmax = ratpmax;


    % Y min P/ATP bei ratpmin
    fixedFluxes(idx_ratp) = ratpmin; % ATP Maintanance
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, product, atpv,fixedFluxes,[]);
    prody(i).YmaxPA_ramin = maxyield;
    prody(i).sup(7) = flux_vec(idx_rsub);
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, -product, atpv,fixedFluxes,[]);
    prody(i).YminPA_ramin = -maxyield;
    prody(i).sup(8) = flux_vec(idx_rsub);
    
    % Y min P/S bei ratpmin
    fixedFluxes(idx_ratp) = ratpmin;
    [maxyield,flux_vec1,success, ~]= CNAoptimizeYield(cnap, -product, substrate,fixedFluxes,[]);
    prody(i).YminPS_ramin = -maxyield;
    prody(i).sup(9) = flux_vec(idx_rsub);

    % Y min P/ATP bei ratpmax
    fixedFluxes(idx_ratp) = ratpmax;
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, -product, atpv,fixedFluxes,[]);
    prody(i).YminPA_ramax = -maxyield;
    prody(i).sup(10) = flux_vec(idx_rsub);
    
    % Y min P/S bei ratpmax
    fixedFluxes(idx_ratp) = ratpmax;
    [maxyield,flux_vec3,success, ~]= CNAoptimizeYield(cnap, -product, substrate,fixedFluxes,[]);
    prody(i).YminPS_ramax = -maxyield;
    prody(i).sup(11) = flux_vec(idx_rsub);
    
        
    % Y min P/A
    fixedFluxes(idx_ratp) = NaN;
    fixedFluxes(idx_rsub) = NaN;
    [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, -product, atpv,fixedFluxes,[]);
    prody(i).YminPA = -maxyield;
    prody(i).ratpPAmin = flux_vec(idx_ratp);
    prody(i).sup(12) = flux_vec(idx_rsub);
    % Y min P/S bei Y min P/A
    fixedFluxes(idx_ratp) = flux_vec(idx_ratp);
    fixedFluxes(idx_rsub) = flux_vec(idx_rsub);
    [maxyield,flux_vec2,success, ~]= CNAoptimizeYield(cnap, -product, substrate,fixedFluxes,[]);
    prody(i).YminPSratpPAmin = -maxyield;
    
    fixedFluxes(idx_ratp) = NaN;
    fixedFluxes(idx_rsub) = NaN;
    
    % Y min P/S, ratp = 0
    fixedFluxes(idx_ratp) = 0;
    [maxyield,flux_vec0,success, ~]= CNAoptimizeYield(cnap, -product, substrate,fixedFluxes,[]);
    prody(i).YminPS_ra0 = -maxyield;
    prody(i).sup(13) = flux_vec(idx_rsub);

    %[yieldspace,success] = CNAplotYieldSpace(cnap,muev,substrate,product,substrate,20,fixedFluxes,[]);

    % plot yieldspace (no longer needed)
    % fill([yieldspace(1:2:end,2); flip(yieldspace(2:2:end,2))],[yieldspace(1:2:end,1); flip(yieldspace(2:2:end,1))],'black');

%     % Y max P/S
%     [maxyield,flux_vec,success, ~]= CNAoptimizeYield(cnap, product, substrate,fixedFluxes,[]);
%     prody(i).sup(1) = flux_vec(idx_rsub);
%     prody(i).YmaxPS = maxyield;
    


end
    clear idx_sprod map_comp_orig idx_wtbm KOreactions reactions...
        KOreactions1 cnap2 fixedFluxes minmax product output ...
        idx_rprod idx_rsub substrate idx_KOreacs map_comp_orig  mue i...
        ratpmax ratpmin prod pname optval muev muer mcs_orig_reacID maxyield...
        idx_ratp idx_rmue flux_vec flux_vec0 flux_vec1 flux_vec2 flux_vec3 
    % output_short mcs
prody_short = prody(find(~cellfun(@isempty,{prody.YmaxPS}')));
prody_short = prody_short(~isnan([prody_short.YminPA_ramin]));
plotyields_one_task;
cnap = cnapbu;