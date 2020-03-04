function [ gMCS , rMCS, sMCS ] = CNAgetGeneMCS( cnap, mcs )
% function 'getGeneMCS'
% ---------------------------
% Trace back how a reaction cut set could be established by gene cuts.
% Returns all possible options plus the resulting 'bigger' reaction cut
% sets, if other reactions are affected by the gene cuts.
% 
% Usage: [ gMCS , rMCS, sMCS ] = getGeneMCS( cnap, mcs )
%
% Input: 
%           cnap:   CellNetAnalyzer project. (Reaction must have a field
%                   named 'geneProductAssociation' in the reaction notes that represent
%                   the gene-reaction rules.)
%           mcs:    [n x r double] matrix that contains MCS
%                   n: number of MCS, r: number of reactions
%
% Output:   
%           gMCS:   {n cell {z cell {g [char]}}} contains all corresponding 
%                   gene cut sets. A nested cell array gMCS{n}{z}{g} with 
%                   n: number of MCS, z: number of possible gene-cut-sets 
%                   that lead to the initial MCS, g: number of genes in the
%                   gene-cut-set
%           
%           rMCS:   {n cell [z x r double]} contains reaction cut sets with
%                   additional cuts that are necessary when gene cuts are 
%                   applied.
%
%           sMCS:   (n struct) contains all information together. Each
%                   position in the struct corresponds to one original MCS.
%                   It lists all possible gene-cut-sets and their extended 
%                   reaction cut sets and also if the cut-set is identical
%                   to the original MCS

% How does it work:
% 1. Extracts all gene-reaction rules
% 2. Find all genes that are associated with reaction cuts in any and all MCS
% 3. Create a small gene-reaction network that contains all associations.
% 4. Conjoin the reaction rules of the reactions for a certain cut set and 
%    find all minimum combinations that lead to the disruption of all reactions 
%    in the cut set. i.e. 3 reactions:     (a | b) | (c & d) | e !=0
%                         <=>   ~a & ~b & ~c & ~e | ~a & ~b & ~d & ~e == 1
% 5. Simulate the gene-kock-outs in the gene-reaction-network and see what
%    is the minimum set of reactions that is disabled through the gene-KOs

%% 1. Extract gene-reaction rules
    geneRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    %geneRules = strrep(strrep(geneRules,'and','&'),'or','|');

    [cnap,~,genes] = CNAgenerateGERassociation(cnap); % get a list of all genes

    % for all reactions, list genes that occur in the gene rule
    affiliatedGenes = cellfun(@(x) unique(strsplit(x)),strrep(strrep(geneRules,'(',''),')',''), 'UniformOutput' , false);
    affiliatedGenes = cellfun(@(x) x(ismember(char(x{:}),genes)),affiliatedGenes,'UniformOutput', false);

%% 2. Find all genes that could occur and find all reactions that could possibly be affected by these genes
    occReacCuts = any(mcs,1);
    occGenes    = unique([affiliatedGenes{occReacCuts}]);
    affectReac  = find(cell2mat(cellfun(@(x) any(ismember(x,occGenes)),affiliatedGenes,'UniformOutput',false)));

%% 3. build a reduced gene-reaction model
    A = eye(length(occGenes));      % Matrix
    m = strcat('spec_',occGenes');  % species
    q = occGenes';                  % reactions
    type = repmat('g',1,length(occGenes)); % Type of reaction: source term genes, sink term reaction, anyth. inbetween
    enzReac = {};
    % add gene-pseudo-enzyme-relations, therefore get gene rule and model the rule as
    % reaction in the network
    for i = affectReac'
        if contains(char(geneRules(i)),'or') % necessary, because the 'children' function would split at '&' otherwise
            gR = children(evalin(symengine,char(geneRules(i))));
        else
            gR = evalin(symengine,char(geneRules(i)));
        end
        gR = arrayfun(@(x) strsplit(strrep(char(x),'&','and'),' and '),gR,'UniformOutput' , false);
        % Only add rule if none of the gene rules are independent from the
        % possible cuts. If there is one rule that relies only on 'other'
        % genes, the reaction would never be really cut through a gene cut
        % set of the original mcs.
        isIndependentGeneRule = cell2mat(cellfun(@(x) ~any(ismember(x,occGenes)),gR,'UniformOutput' , false));
        if ~any(isIndependentGeneRule)
            enzReac(affectReac==i) = {size(A,1)+1:size(A,1)+sum(~isIndependentGeneRule)};
            % for each reaction this loop is run n times. n ist the number
            % of disjuncted ('|') terms in the reaction rule.
            for j = 1:length(gR)
                % Don't add all gene from the rule, only the important
                % ones. The genes that could be cut to genenerate the MCS
                idx = cellfun(@(x) findStrPos(occGenes,x),gR{j},'UniformOutput' , false); 
                idx = [idx{:}];
                % make the gene-peudo-enzyme rule. One column represents a
                % conjuncted ('&') subterm.
                A(:,end+1)      = -iv(size(A,1),idx);
                A(end+1,end)    = 1;
                % name and type
                m = [m; strcat('spec_', cellstr(cnap.reacID(i,:)),'_enz_',num2str(j))];
                q = [q; strcat(cellstr(cnap.reacID(i,:)),'_enz_',num2str(j))];
                type(end+1) = 'a';
            end
        end
    end
    % add pseudo-enzyme-reaction-association
    for i = affectReac'
        if ~isempty(enzReac{affectReac==i})
            m = [m; strcat('spec_', cellstr(cnap.reacID(i,:)))];
            % if there are '&' operators in the gene rule, multiple
            % pseudo-enzymes are associated to one reaction (so multiple
            % columns are needed in the association matrix)
            for j = enzReac{affectReac==i}
                A(:,end+1)      = -iv(size(A,1),j);
                A(size(m,1),end)= 1;
                q = [q; strcat(cellstr(cnap.reacID(i,:)),'_reac_',num2str(j))];
                type(end+1) = 'e';
            end
        end
    end
    % add reaction sinks
    reacs = ~cellfun(@isempty,enzReac);
    reacs = affectReac(reacs);
    numr  = length(reacs);
    A(end-numr+1:end,end+1:end+numr) = -eye(numr);
    q = [q; cellstr(cnap.reacID(reacs,:))];
    type(end+1:end+numr) = 'r';
    %% Build CPlex problem
    cpx = Cplex();
    cpx.DisplayFunc = [];
    cpx.Model.A = A;
    cpx.Model.lb  = zeros(1,size(A,2));
    cpx.Model.ub  =  ones(1,size(A,2));
    cpx.Model.ub(type=='g') =  sum(A(type == 'g',:)<0,2); % upper limit: number of gene Rules in which the gene occurs
    cpx.Model.obj =  double(type=='r'); % maximize number of available reactions
    cpx.Model.sense = 'maximize';
    cpx.Model.ctype(type=='r') = 'B';
    cpx.Model.ctype(type~='r') = 'I';
    cpx.Model.lhs = zeros(1,size(A,1));
    cpx.Model.rhs = zeros(1,size(A,1)); % sollte das nicht auch 0 sein? -> geändert. War vorher Inf
    cpx.solve;
    sol = cpx.Solution.x(end-numr+1:end); % check only reactions
    if ~all(sol)
        error('something''s wrong with the gene-reaction association or the CPLEX netwok');
    end

    %% 4. Investigate all MCS and see whether their
    rMCS = cell(size(mcs,1),1);
    gMCS = cell(size(mcs,1),1);
    for i = 1:size(mcs,1)
        % prepare (find reaction cuts with gene association, i.e. NOT EX_o2_e)
        CutsGene      = find(~cellfun(@isempty,affiliatedGenes)&mcs(i,:)');  % reactions with gene rules
        CutsNoGene    = find(cellfun(@isempty,affiliatedGenes)&mcs(i,:)');
        % generate gene-rule that represents the given cut-set
        geneRule = strjoin(strcat('(',geneRules(CutsGene),')'),' or ');
        % expand rule to DNF
        if contains(char(geneRule),'and') && contains(char(geneRule),'or')
            subterms = children(expand(~evalin(symengine, geneRule)));
        else
            subterms = ~evalin(symengine,geneRule);
        end
        % extract the genes from the expression
        gMCS{i} = arrayfun(@(x) strsplit(strrep(strrep(strrep(char(x),'~','not '),'&','and'),'not ',''),' and '),subterms,'UniformOutput' ,false);
        newRMCS = zeros(size(gMCS{i},2),cnap.numr);
        % simulate the gene knock out and check which reactions can no
        % longer be active
        for j = 1:size(gMCS{i},2)
            cpx.Model.ub(type=='g') =  sum(A(type == 'g',:)<0,2); % upper limit: number of gene Rules in which the gene occurs
            cpx.Model.ub(findStrPos(q,gMCS{i}{j})) = 0;
            cpx.solve;
            sol = cpx.Solution.x(end-numr+1:end); % check only reactions
            newRMCS(j,1:cnap.numr) = iv(cnap.numr,[reacs(~sol);CutsNoGene])';
            if any(mcs(i,:)>iv(cnap.numr,[reacs(~sol);CutsNoGene])')
                disp('kann nicht sein');
            end
        end
        %newRMCS = unique(newRMCS,'rows');
        rMCS(i) = {newRMCS};
    end
    % remove duplicates and restructure stuff
    for i = 1:size(rMCS,1)
        urMCS = unique(rMCS{i},'rows','stable');
        sMCS(i).mcs = mcs(i,:);
        for j = 1:size(urMCS)
            sMCS(i).rMCS(j).rMCS = urMCS(j,:);
            uniqidx = ismember(rMCS{i},urMCS(j,:),'rows');
            sMCS(i).rMCS(j).gMCS = gMCS{i}(uniqidx);
            sMCS(i).rMCS(j).identWithOrig = all(sMCS(i).rMCS(j).rMCS == mcs(i,:));
        end
    end

end

