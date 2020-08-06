[cnap,enzymes] = CNAgenerateGERassociation(iJOcore);
grRules = CNAgetGenericReactionData_as_array(iJOcore,'geneProductAssociation');
%% compress first
[~,~,map11,~,cmp_cnap]=CNAcompressMFNetwork(iJOcore,[],[],1,0,1,[],1);
% restructure gene rules
newrule = repmat({''},1,cmp_cnap.numr);
for i = 1:cmp_cnap.numr
    rules = grRules(logical(map11(:,i)));
    rules = rules(~cellfun(@isempty,rules));
    newrule{i} = strjoin(rules, ') and (');
    if length(rules) > 1
        newrule{i} = ['(' newrule{i} ')'];
    end
end
newrule = strrep(newrule, 'and', '&');
newrule = strrep(newrule, 'or', '|');
syms(cnap.genes)
for i = 1:cmp_cnap.numr
    if ~isempty(newrule{i})
    	newrule{i} = char(expand(str2sym(newrule{i}))); % expand gene rule
    end
end
cellfun(@clear,cnap.genes);
clear rules i;
newrule = strrep(newrule, '&', 'and');
newrule = strrep(newrule, '|', 'or' );
cmp_cnap = CNAsetGenericReactionData_with_array(cmp_cnap,'geneProductAssociation',newrule);
% extend model
cnap_ger = CNAextendModelGenesEnzymes(cmp_cnap);
knockable = cnap_ger.rGeEnReIndicVec == 'g';
% compress again
[~,~,map12,~,cmp_cnap_ger_1]=CNAcompressMFNetwork(cnap_ger,[],[],1,0,1,[],1);

%% append first
cnap_ger = CNAextendModelGenesEnzymes(iJOcore);
knockable = find(cnap_ger.rGeEnReIndicVec == 'g');
[~,~,map21,~,cmp_cnap_ger_21]  = CNAcompressMFNetwork(cnap_ger,knockable,[],1,1,1,[],1);
[~,~,map22,~,cmp_cnap_ger_22] = CNAcompressMFNetwork(cmp_cnap_ger_21,[],[],1,0,1,[],1);