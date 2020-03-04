function [cmp_cnap,gkoCost,gkiCost] = compressGeneRules(cmp_cnap,grRules_old,genes,compMap,gkoCost,gkiCost)
grRules_old = grRules_old(:);
%% translate rules and find identical and contained rules
    rule = cell(length(grRules_old),1);
    for i = 1:length(grRules_old)
        textRule = strrep(strrep(grRules_old{i},'(',''),')','');
        % dissect summands
        disjuctTerms = strsplit(textRule,'or');
        for j = 1:length(disjuctTerms)
            % dissect factors
            factors = strtrim(strsplit(disjuctTerms{j},'and'));
            rule{i}(:,j) = ismember(genes,factors);
        end
    end
    emptyRules = cellfun(@isempty,grRules_old);
    ruleMat{i} = repmat({[]},1,cmp_cnap.numr);
    for i = 1:cmp_cnap.numr
        reacs = find(logical(compMap(:,i)) & ~emptyRules); % find reaction rules to join
        if ~isempty(reacs)
            ruleMat{i} = rule{reacs(1)}; % start with reaction rule of first of the lumped reactions
            for j = 2:length(reacs) % multiply rule for every AND that occurrs (equivalent to expansion of (a+b)*c = a*c + b*c)
                ruleMat{i} = repelem(ruleMat{i},1,size(rule{reacs(j)},2)) | repmat(rule{reacs(j)},1,size(ruleMat{i},2));
                if size(rule{reacs(j)},2) > 1 || j == length(reacs) % compress rules everytime they multiply
                    ruleMat{i} = unique(ruleMat{i}','rows')';
                    ismin = false(1,size(ruleMat{i},2));
                    for k = 1:size(ruleMat{i},2) % mark and eliminate redundand and non-minimal gene rules:
                        ismin(k) = sum(all(ruleMat{i} >= ruleMat{i}(:,k)))==1; % (A and B) or B -> B
                    end                                                        % thus rule (A and B) can be eliminated
                    ruleMat{i} = ruleMat{i}(:,ismin); % only the 'minimal' rules are kept (B)
                end
            end
        end
    end
    str = repmat({{''}},cmp_cnap.numr,1);
    for i = 1:cmp_cnap.numr
        for j = 1:size(ruleMat{i},2)
            str{i}{j} = strjoin(genes(ruleMat{i}(:,j)),' and ');
            if size(ruleMat{i},2) > 1 && sum(ruleMat{i}(:,j)) > 1
                str{i}{j} = ['(' str{i}{j} ')'];
            end
        end
        grRules_new{i} = strjoin(str{i},' or ');
    end
    cmp_cnap = CNAsetGenericReactionData_with_array(cmp_cnap,'geneProductAssociation',grRules_new);
end