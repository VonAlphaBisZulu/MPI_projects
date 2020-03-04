function rmcs = gmcs2rmcs(gmcs,enzymes,rmap,rType)
%
% function 'gmcs2rmcs'
% ---------------------------------------------
% --> translates gene cut sets to reaction cut sets
%
% Usage: rmcs = gmcs2rmcs(gmcs,enzymes,rmap,rType)
%
%  gmcs: provide as: num_genes x num_mcs with following values as intervention indicators: 
%                    cuts: -1, gene additions: 1, not-added genes: nan, else: 0
%
%  enzymes: describes the gene-enzyme-reaction association of the network; it is a
%      struct vector where the i-th entry enzymes(i) represents the i-th enzyme; the fields
%      'genes' and 'reactions' of each entry contain vectors with indices of the genes and reactions
%      associated with this enzyme i; all of these fields must be non-empty, i.e. each enzyme
%      must catalyze at least one reaction and must be associated with at least one gene
%
%  rmap: a vector that maps the reactions from the genome-included cnap to the original 
%
% The following results are returned:
%
%  rmcs: provided as: num_reacs x num_mcs
% 
% Philipp Schneider
% 17.07.2019
%% translate structure enzymes-array (n-1-m) to gr_rules-array (1-1-m) if necessary
if isfield(enzymes,'reactions')
    gr_rules = repelem(enzymes,cellfun(@length,{enzymes.reactions}));
    count = repelem(1:length(enzymes),cellfun(@length,{enzymes.reactions}));
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(count)==x).*(abs(count)==x),unique(count),'UniformOutput',0)'),1)';
    for i = 1:length(gr_rules)
        gr_rules(i).reaction = gr_rules(i).reactions(occ_count(i));
    end
    gr_rules = rmfield(gr_rules,'reactions');
else
    gr_rules = enzymes;
end

%% translate reac MCS
if any(rType=='r')
    ext2old = arrayfun(@(x) find(rmap(x,:),1),1:size(rmap,1));
    rmcs2 = gmcs(ext2old,:);
end

%% translate geneMCS
numr = size(rmap,1);
gmcs = gmcs(rType == 'g',:);

% identified knocked-out or knocked-in rules
rmcs1 = zeros(numr,size(gmcs,2));
for i = 1:size(gmcs,2)
    rule_ki_ko = zeros(1,length(gr_rules));
    for j = 1:length(gr_rules)
        if any(gmcs(gr_rules(j).genes,i) == -1) && ~any(gmcs(gr_rules(j).genes,i) == 1) && all(~isnan(gmcs(gr_rules(j).genes,i)))
            rule_ki_ko(j) = -1;
        elseif all(gmcs(gr_rules(j).genes,i) == 1) && ~any(isnan(gmcs(gr_rules(j).genes,i)))
            rule_ki_ko(j) = 1;
        elseif any(isnan(gmcs(gr_rules(j).genes,i)))
            rule_ki_ko(j) = nan;
        end
    end
    % trace back knocked-out or knocked-in reactions
    rlist = [gr_rules(:).reaction];
    for j = unique(rlist)
        ko_ki_state = rule_ki_ko(j == rlist);
        if all(ko_ki_state == -1)
            rmcs1(j,i) = -1;
        elseif any(ko_ki_state == 1)
            rmcs1(j,i) = 1;
        elseif any(isnan(ko_ki_state)) && ~any(ko_ki_state == 1)
            rmcs1(j,i) = nan;
        elseif all(ko_ki_state == 0 | ko_ki_state == -1)
            rmcs1(j,i) = 0;
        else
            warning('this should not be reached');
        end
    end
end
if any(rType=='r')
    rmcs = rmcs1+rmcs2;
else
    rmcs = rmcs1;
end
end
