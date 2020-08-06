if ~exist('cnan','var') 
    startcna(1); 
end

if ~exist('stored','var') 
    load('_StrainBooster/_My_Simulations/Solutions/substrates.mat','stored');
end
if ~exist('substreacs','var') 
    load('_StrainBooster/_My_Simulations/Solutions/substrates.mat','substreacs');
end

%% Load SBML-Model if necessary
% idx_prod = strcmp(cellstr(iJOBU.reacID),'R_EX_cit_LPAREN_e_RPAREN_');

% if ~exist('iJOBU','var') 
%     iJOBU = TranslateSBML('_StrainBooster/_My_Simulations/iJO1366.xml');
%     [iJOBU, errval] = CNAsbmlModel2MFNetwork(iJOBU);
%     model = 'iJO1366';
%     if errval
%         error('sbml model could not be imported');
%     end
% end
    
k = {};
%% Iterate through product list / models that were used for MCS-calculation
for iProd = 1:length(stored)
    iJOBU = stored(iProd).iJO;
    idx_prod = strcmp(cellstr(iJOBU.reacID),['R_EX_' stored(iProd).name '_LPAREN_e_RPAREN_']);
    
    % boundaries
    iJOBU.reacMax([7 1215]) = 0; % Biomass WT,FHL
    iJOBU.reacMax(9:332) = 0; % Sinks
    iJOBU.reacMax([36,74,82,85,86,95,124,127,128,138,185,186,...
        187,206,208,228,233,237,239,244,245,252,263,291,293,332]) = Inf;
    iJOBU.reacMin(7) = 0; % Biomass WT
    iJOBU.reacMin(1859) = -Inf; % NADH16pp
    iJOBU.reacMin(9:332) = 0; % Sources
    iJOBU.reacMin([74,82,85,86,95,127,128,185,187,206,233,...
        237,239,244,245,252,263,291,332]) = -Inf;
    iJOBU.reacMax(iJOBU.reacMax== 1000) =  Inf;
    iJOBU.reacMin(iJOBU.reacMin==-1000) = -Inf;
    iJOBU.epsilon = 1e-7;
    
    substy = [];
    %% get maximum theoretical yield for all substrates
    for i = 1:size(substreacs,1)

        iJO = iJOBU;
        idx_subs(i) = find(strcmp(cellstr(iJO.reacID),char(substreacs(i,2))));

        substrate               = zeros(1,size(iJO.reacID,1));
        substrate(idx_subs(i))  = -1;
        iJO.reacMin(idx_subs(i))= -10;
        iJO.reacMax(idx_subs(i))= 0;
        
        
        product             = zeros(1,size(iJO.reacID,1));
        product(idx_prod)   = 1; % positive = yield maximization, negative = minimization
        iJO.reacMin(idx_prod) = 0;
        iJO.reacMax(idx_prod) = Inf;

        [maxyield,flux_vec1,success1, ~]= CNAoptimizeYield(iJO, product, substrate,[],[]);
        substy(i) = maxyield;
    end

    %% Find Optimum medium composition. Prefer Glucose only

    iJO = iJOBU;

    iJO.objFunc(:) = 0;

    for i = 1:length(idx_subs)
        iJO.reacMax(idx_subs(i))= 0;
        if isnan(substy(i))
            continue;
        end
        iJO.reacMin(idx_subs(i))= -10;
        switch idx_subs(i)
            case 164 % Glucose (Give Glucose a slight advantage over the other substrates)
                iJO.objFunc(idx_subs(1)) = -substy(1)+iJOBU.epsilon;
            otherwise % Other Substrates
                iJO.objFunc(idx_subs(i)) = -substy(i)-iJOBU.epsilon;
        end
    end

    iJO.reacMin(idx_prod) = 0;
    iJO.reacMax(idx_prod) = Inf;
    iJO.objFunc(idx_prod) = -1;

    [flux_vec2, success2, ~, ~]= CNAoptimizeFlux(iJO);

    %% Output

    str=[];
    idx_rmue = 8;
    for reacIndex = 1:iJO.numr
        educts='';
        zw=find(iJO.stoichMat(:,reacIndex)<0);

        if~isempty(zw)
            educts=[num2str(-iJO.stoichMat(zw(1),reacIndex)),' ',deblank(iJO.specID(zw(1),:))];
            for j=2:length(zw)
                educts=[educts,' + ',num2str(-iJO.stoichMat(zw(j),reacIndex)),' ',deblank(iJO.specID(zw(j),:))];
            end
        end


        products='';
        zw=find(iJO.stoichMat(:,reacIndex)>0);
        if~isempty(zw)
            products=[products,num2str(iJO.stoichMat(zw(1),reacIndex)),' ',deblank(iJO.specID(zw(1),:))];
            for j=2:length(zw)
                products=[products,' + ',num2str(iJO.stoichMat(zw(j),reacIndex)),' ',deblank(iJO.specID(zw(j),:))];
            end
        end

        str=[str; cellstr([educts ' => ' products])];
    end

    % output       = [cellstr(iJO.reacID) num2cell(flux_vec1) str];
    % output_short = output(cell2mat(output(:,2))~=0 & ~isnan(cell2mat(output(:,2))),:);
    % output_very_short = output(abs(cell2mat(output(:,2)))>0.05  & ~isnan(cell2mat(output(:,2))),:);

    % disp(output_very_short);
    % disp(maxyield);

    output2       = [cellstr(iJO.reacID) num2cell(flux_vec2) str];
    output_short2 = output2(cell2mat(output2(:,2))~=0 & ~isnan(cell2mat(output2(:,2))),:);
    output_very_short2 = output2(abs(cell2mat(output2(:,2)))>0.05  & ~isnan(cell2mat(output2(:,2))),:);

    % disp(output_very_short2);
    % Show active substrates
    substFound = idx_subs(ismember(idx_subs,find(cell2mat(output2(:,2))<-1e-5)));
    k(1:length(substFound),iProd) = cellstr(iJO.reacID(substFound,:));
end

% minmax = 'min';
% if sum(product) > 0
%     minmax = 'max';
% end
% disp([minmax ' yield = ' num2str(maxyield) ...
%     ' (' char(strtrim(iJO.reacID(find(product),:))) ' / ' strtrim(char(iJO.reacID(find(substrate,1,'first'),:))) ')' char(10)...
%     num2str(size(output_short,1)) ' active reactions.' char(10)...
%     'mue = ' num2str(output{idx_rmue,2}) char(10)...
%     'product yield = ' num2str((2*(sum(product)<0)-1)*cell2mat(output(idx_prod,2))/cell2mat(output(substidx(1),2)))...
%     ' (' strtrim(char(iJO.reacID(find(idx_prod),:))) ' / ' strtrim(char(iJO.reacID(find(substrate,1,'first'),:))) ')' char(10)...
%     strtrim(char(iJO.reacID(find(idx_prod),:))) ' = '  num2str(cell2mat(output(idx_prod,2)))]);