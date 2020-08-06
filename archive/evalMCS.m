if ~exist('cnan','var')
    startcna(1);
end

if exist('cprodidx','var')
    if isnumeric(cprodidx)
        cprodidx = num2str(cprodidx,'%02i');
    end
    if isstr(cprodidx)
        if cprodidx(1) == 'b'
            cprodidxStr = ['.*2StepCalc.*P(' cprodidx(2:end) ')_.*(?<!cMCS)(?=\.mat)'];
        else
            cprodidxStr = ['.*(?<!2StepCalc-)P(' cprodidx ')_.*(?<!cMCS)(?=\.mat)'];
        end
    end
else
    error('Specify product in variable cprodidx');
end

% find solution files in path
path = './_StrainBooster/_My_Simulations/Solutions/';
filesinPath = dir(path);
%prod = regexp({filesinPath.name}, '.*P([0-6][0-9])_.*(?<!cMCS)(?=\.mat)', 'match'); % PXX_ABC
%prod = regexp({filesinPath.name}, ['.*P(' cprodidx ')_.*(?<!cMCS)(?=\.mat)'], 'match');
prod = regexp({filesinPath.name}, cprodidxStr, 'match');
prod = prod(~cellfun(@isempty,prod));

% prepare compressed iJO
cnapComp = CNAloadNetwork({'iJOcore';1},1,1);

for i = 1:length(prod)
    load([path char(prod{i}) '.mat'],'mcs','cnap','T','D');
    idx.prod = find(xor(T(1,:)',strcmp(cellstr(cnap.reacID),'EX_glc__D_e')));
    idx.subs = find(strcmp(cellstr(cnap.reacID),'EX_glc__D_e'));
    idx.bm   = find(D(1,:));
    idx.atpm = find(D(2,:));
    idx.o2   = find(strcmp(cellstr(cnap.reacID),'EX_o2_e'));
    
    if isempty(idx.o2)
        idx.o2 = find(strcmp(cellstr(cnap.reacID),'R_EX_o2_LPAREN_e_RPAREN_'));
    end
    if length(idx.prod)==2
        idx.prod = find(xor(T(1,:)',strcmp(cellstr(cnap.reacID),'R_EX_glc_LPAREN_e_RPAREN_')));
        idx.subs = find(strcmp(cellstr(cnap.reacID),'R_EX_glc_LPAREN_e_RPAREN_'));
        cnap.reacMin(idx.subs) = -10;
    end
    
    if ~isempty(mcs)
        if size(mcs,1)>1000
            mcs = mcs(unique(round(linspace(1,size(mcs,1),1000))),:);
            warning('only 1000 MCS will be evaluated. Selected through linspace');
        end
        characterizeCMCSxls(cnap,cnapComp,mcs,idx,[path char(prod{i}) '_core' '.xls']);
    else
        warning('no mcs for product.');
        cell2csv([path char(prod{i}) '_core' '.xls'],{'no mcs'});
    end
end