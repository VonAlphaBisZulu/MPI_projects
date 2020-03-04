%path = './_StrainBooster/_My_Simulations/ottoSolutions/';
path = './_StrainBooster/_My_Simulations/Solutions/';
filesinPath = dir(path);
prod = regexp({filesinPath.name}, '.*P([0-6][0-9])_.*\.mat', 'match'); % PXX_ABC
prod = prod(~cellfun(@isempty,prod));
validateMCS = 0;

for j = 1:length(prod)
    solpath = [path char(prod{j})];
        load(solpath, 'iJO', 'mcs','T');
        iJOBU = iJO;
        mcsValid = [];
        numMCS = size(mcs,1);
        numReacs = size(mcs,2);
        prodIndex = find(xor(T(1,:)',strcmp(cellstr(iJO.reacID),'R_EX_glc__D_e')));
        glcIndex = find(strcmp(cellstr(iJO.reacID),'R_EX_glc__D_e'));
        %disp(iJO.reacID(prodIndex,:));
    if numMCS ~= 0
        if numMCS > 100
            mcsIndices = round(linspace(1,size(mcs,1),50));
        else
            mcsIndices = 1:numMCS;
        end
        mincuts = -1;
        if validateMCS
            for i = mcsIndices
                iJO = iJOBU;
                iJO.reacMin(find(mcs(i,:))) = 0;
                iJO.reacMax(find(mcs(i,:))) = 0;
                yields = [-CNAoptimizeYield(iJO,-iv(numReacs,prodIndex)',-iv(numReacs,glcIndex)') CNAoptimizeYield(iJO,iv(numReacs,prodIndex)',-iv(numReacs,glcIndex)')];
                disp(yields);
                if ~any(isnan(yields))
                    mcsValid(end+1,:) = mcs(i,:);
                end
            end
        mincuts = min(sum(repmat(~strcmp(cellstr(iJO.reacID),'R_EX_o2_e'),1,size(mcsValid,1))'.*mcsValid,2));
        end
        if mincuts ~= min(sum(repmat(~strcmp(cellstr(iJO.reacID),'R_EX_o2_e'),1,size(mcs,1))'.*mcs,2));
            if validateMCS
                disp('the smallest MCS was not in validation pool');
            end
            disp([char(prod{j}) ' min cuts: ' char(10) num2str(min(sum(repmat(~strcmp(cellstr(iJO.reacID),'R_EX_o2_e'),1,size(mcs,1))'.*mcs,2)))]);
        end
        if validateMCS
            save(solpath,'mincuts','mcsValid','-append')
        end
    else
        disp([char(prod{j}) char(10) 'doesnt have any cMCS']);
    end
end