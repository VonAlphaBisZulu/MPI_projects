function [cnap, enzymes, genes, gr_rules] = CNAgenerateGERassociation( cnap, grRules, bminimal )
% generates a struct array that contains the information on
% gene-enzyme-reaction-associations. Gene-information can either be
% provided in the cnap-structure as the value of the field 'geneProductAssociation' (see
% CNAsetGenericReactionData), or provided seperately as a cell array.
% Syntax follows Cobra's manner ('and', 'or', curved brackets)
%
% This fuction analyzes the logical structure of the
% gene-reaction-associations to predict "enzymes" or "pseudo-enzymes".
% To detect enzymes correctly it is necessary that all gene-rules are
% fully expanded to DNF. That means no 'or' operator within backets:
%   ORIGINAL:       ((g01 and g02) or (g01 and g03)) and g04
% YES (expanded):    (g01 and g02 and g04) or (g01 and g03 and g04)
% NO (shortest):      g01 and g04 and (g02 or g03)

%% if not provided seperately as cellstr array, read out gene rules from reacNotes
if nargin == 1
    grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    bminimal = 0;
end

if nargin >= 2
    if isempty(grRules)
        grRules = CNAgetGenericReactionData_as_array(cnap,'geneProductAssociation');
    end
end

% enumerate all occurring genes, generate a vector where all of them occur
% only once.
genesStr = {};
for i = 1:length(grRules)
    genesStr = [genesStr; split(regexprep(grRules(i),'\(|\)',''))]; % regexprep removed brackets
end
genesStr = unique(genesStr);
genesStr = genesStr(~strcmp(genesStr,'and')&~strcmp(genesStr,'or')&~strcmp(genesStr,''));

% find all enzymes. Gene rules are evaluated. Genes that are connected
% through an 'and' are concatinated (with '-' as delimiter) and regarded as one enzyme.
enzymeReacAssoc = {};
enzymeGeneAssoc = {};
enzRules = regexprep(regexprep(regexprep(grRules,'  *and  *','-'),'  *or  *',' '),'\(|\)','');
for i = 1:length(grRules)
     eRule = {split(enzRules(i))}; % get seperate enzymes that catalyze this reaction
     for j=1:length(eRule{:})
        enzymeReacAssoc(i,j) = eRule{1}(j); % enter enzmes into a Matrix 
                                            % (reacs -> list of catalyzing enzymes)
        enzymeGeneAssoc(i,j) = {findStrPos(genesStr,split(eRule{1}{j},'-'))}; % enter
                                              % corresponding gene-enzyme-rules at same
                                              % position in the matrix
     end
end
% generate a vector containing all occurring enzymes only once.
enzymesStr = reshape(enzymeReacAssoc,1,size(enzymeReacAssoc,1)*size(enzymeReacAssoc,2));
enzymesStr = unique(enzymesStr(~cellfun(@isempty,enzymesStr)));
enzymesStr = enzymesStr(~strcmp(enzymesStr,''));
enzymes = struct('name',{},'reactions',{},'genes',{},'strReac',{},'strGene',{});

for i=1:length(enzymesStr)
    enzymes(i).name      = strcat('Enz-', enzymesStr(i));
    
    % find enzyme-reaction association
    [reacID,enzInReac]   = find(strcmp(enzymeReacAssoc,enzymesStr(i)));
    enzymes(i).reactions = reacID';
    % find corresponding gene-enzyme association (on 
    enzymes(i).genes     = enzymeGeneAssoc{reacID(1),enzInReac(1)};
    
    % additionally also save gene-enzyme-reaction associations as text
    enzymes(i).strReac   = cellstr(cnap.reacID(enzymes(i).reactions',:))';
    enzymes(i).strGene   = cellstr(genesStr(enzymes(i).genes))';
end

if bminimal
    genesStr = {};
    for i = 1:length(enzymes)
        name = char(enzymes(i).name);
        genesStr(i,:) = [name {enzymes(i).strGene}];
        enzymes(i).strGene = cellstr(name(5:end));
        enzymes(i).genes = i;
    end
end

genes = genesStr;
if ~isempty(enzymes)
    rulecount = repelem(1:length(enzymes),cellfun(@length,{enzymes(:).reactions}));
    countOcc = occurcount(rulecount);
    for i = 1:length(rulecount)
        gr_rules(i).reaction  = enzymes(rulecount(i)).reactions(countOcc(i));
        gr_rules(i).strReac   = enzymes(rulecount(i)).strReac(  countOcc(i));
        gr_rules(i).genes     = enzymes(rulecount(i)).genes;
        gr_rules(i).strGene   = enzymes(rulecount(i)).strGene;
        gr_rules(i).name    = {['Rule-' enzymesStr{rulecount(i)} '-r' num2str(gr_rules(i).reaction)]};
    end
else
    gr_rules = struct('name',{},'reaction',{},'genes',{},'strReac',{},'strGene',{});
end

cnap.enzymes  = enzymes;
cnap.genes    = genesStr;
cnap.gr_rules = gr_rules;

function indices = findStrPos( str , pattern )
% str       the space that is seached
% pattern   the search keyword or pattern
% convert str type to cell
switch class(str)
    case 'string'
        str = cellstr(strtrim(str));
    case 'char'
        str = cellstr(str);
    case 'cell'
        
    otherwise
        error('input 1 of findStrPos doesn''t correct type');
end
% convert pattern type to cell
switch class(pattern)
    case 'string'
        pattern = char(pattern);
    case 'char'
        pattern = cellstr(pattern);
    case 'cell'

    otherwise
        error('input 2 of findStrPos doesn''t correct type');
end
[rows,cols] = size(pattern);
for k = 1:(rows*cols)
    ind                      = find(strcmp(str, pattern(k)));
    indices(1:length(ind),k) = ind;
end
if (any(size(indices)==0))
    indices = [];
end
end
%% occurrency counter [1 2 2 3 3] -> [1 1 2 1 2]
function occ_count = occurcount(y)
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(y)==x).*(abs(y)==x),unique(y),'UniformOutput',0)'),1)';
end
end