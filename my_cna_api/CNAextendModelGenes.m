function [ gcnap, rmap, gmap, gkoCost, gkiCost, gT, gD, rType, sType, gr_rules ] = ...
    CNAextendModelGenes( cnap, gr_rules, koCost, kiCost, T, D, gkoCost, gkiCost, expand_notknockable )
% function 'CNAextendModelGenes'
% ---------------------------
% Extends a given stoichiometric networks with terms for genes and enzymes. 
% Reversible reactions and reactions that can be catalyzed by multiple
% enzymes are copied.
% 
% Usage: [ cnap, notknockable, reacMapOldNew, rGeEnReIndicVec, sGeEnReIndicVec, T, t, D, d ] = ...
%    CNAextendModelGenesEnzymes( cnap, enzymes, notknockable, flagGeneOrEnzymeKO, T, t, D, d, geneEnzNoKO )
%
% Input: 
%           cnap:   CellNetAnalyzer project. (with gene rules in reaction
%                                             notes)
%   <optional>
%               enzymes:            a struct-array that contains the fiels
%                                   reactions [double] and genes [double]
%                                   that represent the gene-enzyme-reaction
%                                   associations
%               koCost, kiCost:     Vectors that specify the knockable reactions.
%                                   Length must be identical to the number of reactions
%                                   of the network. 'nan' -> reaction is not knockable.
%                                   Knock-Ins 'override' knock-outs.
%                                   If reactions with no gene association are marked
%                                   as knockable (such as substrate or O2 uptake), 
%                                   they will keep their 'knockability'.
%                                   In case of enzyme associated reactions, all underlying
%                                   genes (or enzymes) will inherit the KO/KI costs and
%                                   the reaction itself will become not-knockable.
%                                   In case of conflicts, when two reactions share the
%                                   same gene, but one is knockable and the other is not,
%                                   the gene will be set 'not-knockable' / 'nan'.
%               T, t, D, d:         Target and Desired vectors that should
%                                   be reshaped (not yet implemented
%                                   correctly)
%               noKO_gen_enz:       logical vector of notknockable genes (genes only, in case of enzymes,
%                                   the particular enzyme from the enzyme-struct).
%            expand_notknockable:   Add enzyme and gene-association, even if their reactions
%                                   were marked as notknockable (koCost(i) = nan)
%                                   (default: 1)
%
% Output:   
%           cnap:   CellNetAnalyzer project, extended by
%                     1. source reactions for every gene
%                     2. converstion reactions from genes to reaction-pseudo-metabolites defined by the reaction rules
%                     3. pseudo-metabolite consumption is added to catalyzed reactions
%                     - fields sGePmReIndicVec and rGePmReIndicVec that
%                       show whether species or reactions in CNA project
%                       represent real metabolites, genes or reactions (or
%                       gene source terms, translations or chemical
%                       reactions)
%                     - a list of all genes
%           reacMapOldNew:  A matrix that maps the reactions of the old
%                           stoichiometric matrix with the new one

%% get parameters
% if enzmes don't have names, generate names
% if genes don't have names, generate names
% check if cnap has enzyme
if nargin < 2 || isempty(gr_rules)
    if ~isfield(cnap,'enzymes')
        [cnap,~,~,gr_rules] = CNAgenerateGERassociation(cnap);
    else
        gr_rules = cnap.gr_rules;
    end
end
numRules = length(gr_rules);
if nargin < 3 || isempty(koCost)
    koCost = ones(cnap.numr,1);
end
if nargin < 4 || isempty(kiCost)
    kiCost = ones(cnap.numr,1);
end
if nargin < 5
    D = {};
    T = {};
end
if nargin < 7 || isempty(gkoCost) % number of genes is defined by 
    numGenes = max([gr_rules(:).genes]);
    gkoCost = ones(numGenes,1);
    gkoCost_provided = 0;
else
    numGenes = length(gkoCost);
    gkoCost_provided = 1;  
end
if nargin < 8 || isempty(gkiCost)
    gkiCost = nan(numGenes,1);
    gkiCost_provided = 0;   
else, gkiCost_provided = 1;
end
if nargin < 9
    expand_notknockable = 1;
end
c_macro = cnap.macroDefault;

%% 0 initialize
% if enzmes don't have names, generate names
% if genes don't have names, generate names
% check if cnap has enzyme

koCost = double(koCost(:)'); % reshape as row (double)
kiCost = double(kiCost(:)');
koCost(~isnan(kiCost)) = nan; % knock-ins 'override' knock-outs

gkoCost = double(gkoCost(:)'); % reshape as row (double)
gkiCost = double(gkiCost(:)');
gkoCost(~isnan(gkiCost)) = nan; % knock-ins 'override' knock-outs

ruleMat = repmat({[]},1,cnap.numr);
for i = 1:cnap.numr
    rule_idx = find(ismember([gr_rules(:).reaction],i));
    rule_i = arrayfun(@(x) gr_rules(x).genes,rule_idx,'UniformOutput',0);
    ruleMat{i} = sparse(repelem(1:length(rule_idx),cellfun(@length,rule_i)),cell2mat(rule_i),1,length(rule_idx),numGenes);
end

if ~expand_notknockable
    rule_removelist = zeros(1,length(gr_rules));
    [gene_abund,genes] = hist([gr_rules(:).genes],unique([gr_rules(:).genes]));
    if ~gkoCost_provided
        for i = find(isnan(koCost) & ~cellfun(@isempty,ruleMat)) % for notknock reaction with gene reaction rule:
            % if a reaction is notknockable and there are genes that only occur with this reaction,
            % set those genes to notknockable, because it is never advantageous to knock them out.
            genes_1 = [gr_rules([gr_rules(:).reaction] == i).genes];
            for j = unique(genes_1)
                % if the number of occurrencies of this gene in the rules of one reaction is equal to
                % the total number of occurrencies of this gene in all rules.
                if sum(genes_1 == j) == gene_abund(genes == j)
                    gkoCost(j) = nan;
                end
            end
        end
    end
    if ~isempty(gkoCost)
        for i = 1:length(gr_rules)
            % if all genes of one rule are notknockable, delete all rules for the same reaction
            % because the reaction can never be knocked out
            if all(isnan(gkoCost(gr_rules(i).genes))) && all(isnan(gkiCost(gr_rules(i).genes)))
                rule_removelist([gr_rules(:).reaction] == gr_rules(i).reaction) = 1;
            end
        end
    end
    gr_rules = gr_rules(~rule_removelist);
    numRules = length(gr_rules);
end

if isempty(gr_rules)
    error('Model extension with enzymes and genes not possible. Make sure that at least one enzyme can be knocked out.');
end

% set gene and enzyme names
if ~isfield(gr_rules,'strGene')
    geneNames  = strcat('G',num2str(1:numGenes));
else
    [a,b] = unique([gr_rules(:).genes]);
    gname = [gr_rules(:).strGene]';
    geneNames = arrayfun(@(x) ['gene ' num2str(x)],1:numGenes,'UniformOutput',0);
    geneNames(a) = gname(b);
end
if ~isfield(gr_rules,'name')
    ruleNames = strcat('PM',num2str(1:length(gr_rules)));
else
    ruleNames = [gr_rules(:).name]';
    if any(ismember(ruleNames,geneNames))
        ruleNames = strcat('Rule-',ruleNames);
    end
end

% duplicate reactions that have an enzyme association and are reversible and with gene association.
rules_per_reac = hist([gr_rules(:).reaction],1:cnap.numr);
rev_reac_w_rule = rules_per_reac & (sign(cnap.reacMin) == -1)'; % subset of these reactions that are active in reverse direction
for_reac_w_rule = rules_per_reac & (sign(cnap.reacMax) ==  1)'; % subset of these reactions that are active in forward direction
% which reactions need to be split, create mapping.
r_old2new = num2cell(1:cnap.numr);
for i = find(rev_reac_w_rule | for_reac_w_rule) % indicate which reactions are cloned (+/-)
    r_old2new{i} = [i*ones(1,for_reac_w_rule(i)) -i*ones(1,rev_reac_w_rule(i))]; % reactions
end
r_new2old = [r_old2new{:}];
rules_per_reac_new = rules_per_reac(abs(r_new2old));
rType(1:length(r_new2old)) = 'r'; % reaction Type: enzyme catalyzed reactions
sType(1:cnap.nums) = 'm'; % species Type: metabolites
rmap = full(sparse(abs(r_new2old),1:length(r_new2old),sign(r_new2old)));
 % reactions that stem back from reversible reactions with enzymes (and were split or inverted)
 % 0: reaction was not split, 1: positive sense reaction part, -1: reverse sense reaction part
r_was_rev = rev_reac_w_rule*rmap;

%% 1. extend stoichMat
%
% Take original model + gen-enzyme-reaction-association
% Split enzyme catalyzed reactions if necessary:
%   extend stoichMat with k reaction columns to split up reversible
%                          reactions
%   extend stoichMat with l reaction columns to cover cases where
%                          different enzymes can catalyze the same reaction
%
% duplicate reactions that have an enzyme association and are reversible or 
% with isoenzyme association.

% Model: Stoichmat + bounds + objFunc
gcnap.stoichMat = cnap.stoichMat* rmap; % new Stoichmat
gcnap.reacMin = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMax = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMin(~r_was_rev) = cnap.reacMin(abs(r_new2old(~r_was_rev)));
gcnap.reacMin(~~r_was_rev) = 0;
gcnap.reacMax(r_was_rev >= 0) =  cnap.reacMax(abs(r_new2old(r_was_rev >= 0)));
gcnap.reacMax(r_was_rev <  0) = -cnap.reacMin(abs(r_new2old(r_was_rev <  0)));
gcnap.objFunc   = rmap'*cnap.objFunc;
% reacIDs
% Counter for split reactions (only for naming)
gcnap.reacNotes = cnap.reacNotes(abs(r_new2old));
gcnap.reacID    = cellstr(cnap.reacID(abs(r_new2old),:));
gcnap.reacID(rules_per_reac_new>0,:) = strcat(gcnap.reacID(rules_per_reac_new>0,:),'_gen');
gcnap.reacID(r_was_rev==-1,:) = strcat(gcnap.reacID(r_was_rev==-1,:),'_rev');
gcnap.reacID    = char(gcnap.reacID);
% specIDs
gcnap.specID    = cnap.specID;
gcnap.specNotes = cnap.specNotes;

%% 2. add reaction-pseudometabolites and genes
%   extend stoichMat with m1 rows for reaction-pseudometabolites
%                                          (all have values -1 below N
%                                           and 1 to the bottom right of N)
%   extend stoichMat with m2 rows for genes (all bottom right corner)
%   extend stoichMat with n1 = m1 columns for gene-to-reaction-pseudometabolite reactions
%   extend stoichMat with n2 = m2 columns for gene-producing reactions

% create vector that indicates whether reaction is gene expression, gene rule ('p') or reaction
reac_pseudo_met = unique(find(rules_per_reac)); % 1 pseudometabolite for each reaction
numReacPseudom = length(reac_pseudo_met);       %   (not 2 pseudo m. if reaction was split due to rev.)

% add pseudometabolites and pseudometabolites consumption to reactions
sType(end+1:end+numReacPseudom) = 'p';
pseudomet_consump_x = find(ismember(abs(r_new2old),reac_pseudo_met));
pseudomet_consump_y  = repelem(1:numReacPseudom,(rev_reac_w_rule(reac_pseudo_met) & for_reac_w_rule(reac_pseudo_met))+1);
gcnap.stoichMat(sType=='p',rType=='r') = -sparse(pseudomet_consump_y,pseudomet_consump_x,1,numReacPseudom,size(gcnap.stoichMat,2));

% add pseudometabolite source terms (gene consumption is added later)
rType(end+1:end+numRules) = 'p';
pseudomet_source_x = 1:length(gr_rules);
pseudomet_source_y = arrayfun(@(x) find(reac_pseudo_met == x), [gr_rules(:).reaction]);
gcnap.stoichMat(sType=='p',rType=='p') = sparse(pseudomet_source_y,pseudomet_source_x,1,numReacPseudom,numRules); % add pseudometabolite source reaction
gcnap.reacMin(rType=='p') = 0;
gcnap.reacMax(rType=='p') = 1000;
gcnap.specID = char([cellstr(gcnap.specID) ; ...
                     strcat('psmet-',cellstr(cnap.reacID(reac_pseudo_met,:)))]);
gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('GR-',ruleNames)]);

% add genes, gene consumption / gene rule part to pseudometabolite, also add gene source terms
sType(end+1:end+numGenes) = 'g';
rType(end+1:end+numGenes)   = 'g';
% change pseudoenzyme source terms into gene-pseudoenzyme reactions
pseudomet_gene_x = repelem(1:length(gr_rules),cellfun(@length,{gr_rules.genes}));
pseudomet_gene_y = [gr_rules.genes];
gcnap.stoichMat(sType=='g',rType=='p') = -sparse(pseudomet_gene_y,pseudomet_gene_x,1,numGenes,numRules);
gcnap.stoichMat(sType=='g',rType=='g') = -eye(numGenes); % add gene source pseudo-reaction
gcnap.reacMin(rType=='g') = -1000;
gcnap.reacMax(rType=='g') = 0;
unused_genes = setdiff(1:numGenes,[gr_rules(:).genes]); % deactivate unused genes
gcnap.reacMin(find(rType=='g',1)-1+unused_genes) = 0;
gkoCost(unused_genes) = nan;
gkiCost(unused_genes) = nan;
gcnap.specID = char([cellstr(gcnap.specID) ; geneNames(:)]);
gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('SG-',geneNames(:))]);

gcnap.specNotes(size(gcnap.stoichMat,1)) = {''};
gcnap.reacNotes(size(gcnap.stoichMat,2)) = {''};
gcnap.objFunc(rType=='p' | rType=='g')   = 0;
rmap = [rmap, zeros(cnap.numr,numGenes+numRules)];

%% 3. adapt kiCost, koCost, D, d, T, t
% kiCost and koCost vectors will be derived from reaction-Cost vectors
% if reactions have associated enzymes/genes, those will be marked as knockable.
% Else, the reactions will stay knockable. Conflicts are treated as follows:
% When geneA catalyzes R1 (knockable) and R2 (notknockable), the corresponding geneA 
% will be knockable. If R1 has a lower cost than R2, lower cost will be taken into account.

gmap = [zeros(numGenes,sum(rType ~= 'g')), eye(numGenes)];
if gkoCost_provided
    gkoCost = [nan(1,sum(rType ~= 'g')) , gkoCost];
else
    gkoCost = nan(1,length(rType));
end
if gkiCost_provided
    gkiCost = [nan(1,sum(rType ~= 'g')) , gkiCost];
else
    gkiCost = nan(1,length(rType));
end
cgenes = {gr_rules(:).genes};

for i = 1:size(gcnap.stoichMat,2)
    switch rType(i)
        case 'r' % reactions stay knockable if they were knockable before and don't have a gene-association
            if ~rules_per_reac_new(i)
                gkoCost(i) = koCost(r_new2old(i));
                gkiCost(i) = kiCost(r_new2old(i));
            end
        case 'p'
        case 'g'
            if ~gkoCost_provided || ~gkiCost_provided
                gene_id = find(find(rType=='g') == i);
                rules = cellfun(@(x) ismember(gene_id,x),cgenes);
                if any(~isnan(koCost([gr_rules(rules).reaction]))) && ~gkoCost_provided
                    gkoCost(i) = min(koCost([gr_rules(rules).reaction]));
                end
                if any(~isnan(kiCost([gr_rules(rules).reaction]))) && ~gkiCost_provided
                    gkiCost(i) = min(kiCost([gr_rules(rules).reaction]));
                end
            end
        otherwise
            error(['invalid rType. Should be ''r'' for reactions, ''e'''...
                   'for enzyme rules/sources and ''g'' for gene sources']);
    end
end
if any(~isnan(gkiCost) & ~isnan(gkoCost))
    warning('There are conflicts between KIable and KOable vector. Please check.');
end

if ~gkoCost_provided
    knockable_genes = gkoCost(rType=='g');
    not_knockable_warning = '';
    for i = find(isnan(koCost) & rules_per_reac > 0) % check if notknockable reactions still depend on at least one rule that is notknockable
        rules = find([gr_rules(:).reaction] == i);
        independent = 0;
        for j = rules
            if all(isnan(knockable_genes(gr_rules(j).genes)))
                independent = 1;
            end
        end
        if ~independent
                not_knockable_warning = [not_knockable_warning strtrim(cnap.reacID(i,:)) ...
            ' :RULE: ' num2str(rules) ...
            ' :GENE: ' num2str([gr_rules(rules).genes])  '; ' newline];
        end
    end
    if ~isempty(not_knockable_warning)
        warning(['A reaction that was before not knockable has become knockable through the expansion. '...
            'The affected reactions, enzymes and genes are ' newline not_knockable_warning 'Please check for ' ...
            'genetic overlaps of knockable and notknockable reactions before expanding the model.']);
    end
end
% T and D
if ~isempty(T)
    for i = 1:length(T) % fill up to match reaction vector size
        gT{i}  = T{i}*rmap;
    end
else
    gT = {};
end
if ~isempty(D)
    for i = 1:length(D) % fill up to match reaction vector size
        gD{i}  = D{i}*rmap;
    end
else
    gD = {};
end

gcnap = CNAgenerateMFNetwork(gcnap,1);
%% 4. in case of gui, generate map
if cnap.has_gui
    gcnap.local.errval  = 0;
    gcnap.numr          = size(gcnap.stoichMat,2);
    gcnap.nums          = size(gcnap.stoichMat,1);
    gcnap.mue           = find(r_new2old == cnap.mue');
    gcnap.figs          = cnap.figs;
    gcnap.path          = cnap.path;
    gcnap.color1        = cnap.color1;
    gcnap.color2        = cnap.color2;
    gcnap.color3        = cnap.color3;
    gcnap.color4        = cnap.color4;
    gcnap.maps          = cnap.maps;
    gcnap.nummaps       = cnap.nummaps;
    gcnap.show_flux_format= cnap.show_flux_format;
    gcnap.net_var_name  = cnap.net_var_name;
    gcnap.reacDefault   = cnap.reacDefault(abs(r_new2old)).*sign(r_new2old);
    gcnap.reacDefault(r_was_rev == -1) = double(gcnap.reacDefault(r_was_rev == -1)>=0).*gcnap.reacDefault(r_was_rev == -1);
    gcnap.reacFontSize  = cnap.reacFontSize;
    gcnap.reacBoxWidth  = cnap.reacBoxWidth;
    gcnap.reacBoxHeight = cnap.reacBoxHeight;
    gcnap.specFontSize  = cnap.specFontSize;
    gcnap.specBoxWidth  = cnap.specBoxWidth;
    gcnap.specBoxHeight = cnap.specBoxHeight;
    gcnap.specBoxColor  = cnap.specBoxColor;
    %% extend reac Boxes vector
    gcnap.reacBoxes = cnap.reacBoxes(abs(r_new2old),:);
    gcnap.reacBoxes(r_was_rev == -1, 4) = 0; % reversible part
    gcnap.reacBoxes(r_was_rev == -1, 2) = gcnap.reacBoxes(r_was_rev == -1, 2) + 20;
    gcnap.reacBoxes(r_was_rev == -1, 3) = gcnap.reacBoxes(r_was_rev == -1, 3) + 3;
    disp('drawing new reaction boxes');
    for reacIndex = find(r_was_rev == -1)
        % declaration of handle
        cnan.open_projects.(gcnap.net_var_name).gui.handles= struct;
        gcnap = update_after_change(gcnap);
        currmap = gcnap.reacBoxes(reacIndex,5);
        % is visible/editable
        zw=gcnap.reacBoxes(reacIndex,6);
        fig = gcnap.figs(currmap,:);
        % make box
        zw1=uicontrol('Style', 'edit','Parent',fig(1), 'String', '###',...
            'Units','normalized','HorizontalAlignment','left','BackgroundColor',gcnap.color1,'ForegroundColor',gcnap.textColor,'TooltipString',gcnap.reacID(reacIndex,:));
        set(zw1, 'ButtonDownFcn', {@execute_callback, gcnap.net_var_name,...
            {'check_for_right_click', 'reaceditmask'}, {'reacenr', reacIndex}});
        % save handle
        gcnap.reacBoxes(reacIndex,4)=zw1;
        % adjust "zoom"
        place_box(fig,zw1,...
            gcnap.reacBoxes(reacIndex,2),...
            gcnap.reacBoxes(reacIndex,3),...
            gcnap.reacFontSize(currmap),...
            gcnap.reacBoxWidth(currmap),...
            gcnap.reacBoxHeight(currmap));
        % put default rate in textbox
        if isnan(gcnap.reacDefault(reacIndex))
            set(zw1,'String','#');
        else
            set(zw1,'String',num2str(gcnap.reacDefault(reacIndex)));
        end
        % is visible/editable
        if(zw==2)
            set(zw1,'Style', 'text');
        elseif(zw==3)  %%non-visible
            set(zw1,'Visible','off');
        end
    end
    disp('generating map for gene and pseudometabolite reactions');
    gcnap.reacBoxes(rType ~= 'r',5) = -1;
    gcnap = myCNAgenerateMap(gcnap,sum(rType ~= 'r')<=110);
end

% create gcnap Network
gcnap.epsilon = cnap.epsilon;
gcnap.rType = rType;
gcnap.sType = sType;
end

%% supplementary
function place_box(fig,handle,xp,yp,fontsize,box_width,box_height) % copied from zoom_single_box
    zzz=get(fig(1),'CurrentAxes');
    if(~fig(4))
    xxx=get(zzz,'XLIM');
    yyy=get(zzz,'YLIM');
    xxxl=xxx(2)-xxx(1);
    yyyl=yyy(2)-yyy(1);
    xp=(xp-xxx(1))/xxxl;
    yp=(yyy(2)-yp)/yyyl;
    box_width=box_width*fig(3)/xxxl;
    box_height=box_height*fig(2)/yyyl;
    fontsize=fontsize*fig(3)/xxxl;
    end
    pos=[xp yp box_width box_height];
    set(handle,'Position',pos,'FontSize',fontsize);
end
%% occurrency counter [1 2 2 3 3] -> [1 1 2 1 2]
function occ_count = occurcount(y)
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(y)==x).*(abs(y)==x),unique(y),'UniformOutput',0)'),1)';
end
