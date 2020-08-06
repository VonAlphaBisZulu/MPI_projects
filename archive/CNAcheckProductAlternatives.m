function growthWithAltProd = CNAcheckProductAlternatives( cnap, mcs, ProductOI )
    cnapBU = cnap;
    if ischar(ProductOI)
        idx_p = find(strcmp(cellstr(cnap.reacID),ProductOI));
    elseif isnum(ProductOI)
        idx_p = ProductOI;
    else
        error('no valid product specified');
    end
    idx_sp = find(cnap.stoichMat(:,idx_p));
    % Find biomass reaction
    mue = regexp(cellstr(cnap.reacID), '.*BIOMASS.*core.*', 'match');
    mue = find(~cellfun(@isempty,mue));
    cnap.objFunc = -iv(cnap.numr,mue);
    % put names in first column
    growthWithAltProd        = cellstr(cnap.specID);
    growthWithAltProd(end+1) = {cnap.specID(idx_sp,:)};
    growthWithAltProd(end+1) = {'NAD(P)H Recovery'};
    progBar = Progressbar(0,'Testing alternative products',cnap.nums+size(mcs,1)*(2+cnap.nums));
    % 1st test in original model
    for j=1:size(cnap.specID,1)
        cnap = cnapBU;
        cnap = CNAaddReactionMFN(cnap,'auxiliary_sink',[strtrim(cnap.specID(j,:)) ' = '],0,Inf,0, NaN, 0, {}, 1, 1, 0, 0);
        [~,~,~,growth] = CNAoptimizeFlux(cnap);
        growthWithAltProd(j,2) = {-growth};
        progBar = progBar.update(j);
    end
    % 2nd test in mcs scenario
    for i=1:size(mcs,1)
        for j=1:size(cnap.specID,1)
            cnap = cnapBU;
            cnap.reacMin(idx_p) = 0;
            cnap.reacMax(idx_p) = 0;
            cnap.reacMin(find(mcs(i,:))) = 0;
            cnap.reacMax(find(mcs(i,:))) = 0;
            % 2nd test in mcs scenario
            [~,~,~,growth] = CNAoptimizeFlux(cnap);
            growthWithAltProd(j,i+2) = {-growth};
            progBar = progBar.update(i*(cnap.nums+2)+j);
        end
        % test with original product (export on)
        cnap = cnapBU;
        cnap.reacMin(find(mcs(i,:))) = 0;
        cnap.reacMax(find(mcs(i,:))) = 0;
        [~,~,~,growth] = CNAoptimizeFlux(cnap);
        growthWithAltProd(end-1,i+2) = {-growth};
        % test with NADH
        cnap = cnapBU;
        cnap.reacMin(find(mcs(i,:))) = 0;
        cnap.reacMax(find(mcs(i,:))) = 0;
        cnap.reacMin(idx_p) = 0;
        cnap.reacMax(idx_p) = 0;
        cnap = CNAaddReactionMFN(cnap,'nadh_recovery','nadh_c + h_c = nad_c + h2_c',0,Inf,0, NaN, 0, {}, 1, 1, 0, 0);
        cnap = CNAaddReactionMFN(cnap,'nadph_recovery','nadph_c + h_c = nadp_c + h2_c',0,Inf,0, NaN, 0, {}, 1, 1, 0, 0);
        [~,~,~,growth] = CNAoptimizeFlux(cnap);
        growthWithAltProd(end,i+2) = {-growth};
    end
    progBar.delete();
end