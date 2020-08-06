function [stch chrg] = verifyReacStoich( cnap, reacnr )
% check if charge and reaction stoichiometry is feasible
    if ischar(reacnr)
        reacnr = find(strcmp(cellstr(cnap.reacID),reacnr),1);
    end
    chrg = 0;
    stch = 0;
    % read out charge of all metabolites
             % [ spec_idx , stoichFactor ]
    educts   = [find(cnap.stoichMat(:,reacnr)<0) cnap.stoichMat(find(cnap.stoichMat(:,reacnr)<0),reacnr)];
    products = [find(cnap.stoichMat(:,reacnr)>0) cnap.stoichMat(find(cnap.stoichMat(:,reacnr)>0),reacnr)];
    % charge
    formula = regexp(cnap.specNotes, '(?<=\[).*(?=\])' , 'match');
    formula = [formula{:}]';
    if any(cellfun(@isempty,formula))
        error('stoichiometry not specified for all species');
    end
    charge = regexp(cnap.specNotes, '(?<=<).*(?=>)' , 'match');
    if any(cellfun(@isempty,charge))
        error('charge not specified for all species');
    end
    charge = cell2mat(cellfun(@str2num,[charge{:}],'UniformOutput',false));
    
    % Make a matrix with all stoichiometries for species
    atoms = regexp(formula, '[A-Z][a-z]*' , 'match');
    atoms = unique([atoms{:}]);
    stoichTable = zeros(cnap.nums,length(atoms));
    for i = 1:length(atoms)
        atmCount = regexp(formula, [char(atoms(i)) '(?![a-z])[0-9]*'] , 'match');
        atmCount2 = find(~cellfun(@isempty,atmCount));
        atm = repmat({''},cnap.nums,1);
        atm(atmCount2) = [atmCount{atmCount2}];
        one = cellfun(@length,atm)==length(char(atoms(i)));
        atm = cellfun(@(x,i) x(i:end),atm,repmat({length(char(atoms(i)))+1},cnap.nums,1),'UniformOutput',false);
        atm(one) = {'1'};
        atm(cellfun(@isempty,atm)) = {'0'};
        stoichTable(:,i) = cell2mat(cellfun(@str2num,atm,'UniformOutput',false));
    end
    
    % balance atoms
    stoichsum = sum([educts(:,2).*stoichTable(educts(:,1),:);products(:,2).*stoichTable(products(:,1),:)],1);
    stch = [atoms; num2cell(stoichsum)]';
    
    if ~all(cnap.reacID(reacnr,1:3) ~= 'EX_') && any(stoichsum)
        warning('non-sink/source reaction with uneven stoichiometry');
    end
    
    % balance charge
    chrg = sum([educts(:,2).*charge(educts(:,1))';products(:,2).*charge(products(:,1))']);
    end

