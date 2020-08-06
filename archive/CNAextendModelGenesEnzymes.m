function [ gcnap, reacMap_old_new, gkoCost, gkiCost, gT, gD, rType, sType, enzymes ] = ...
    CNAextendModelGenesEnzymes( cnap, enzymes, flagGeneOrEnzymeKO, koCost, kiCost, T, D, gkoCost, gkiCost, expand_notknockable )
% function 'CNAextendModelGenesEnzymes'
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
%               flagGeneOrEnzymeKO: 'g': Genes, 'e': Enzymes (default 'g')
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
%                     2. converstion reactions from genes to pseudo-enzymes defined by the reaction rules
%                     3. enzyme consumption is added to catalyzed reactions
%                     - 'enyzmes'-struct-vector that describes every enzyme's reactions and genes   
%                     - fields sGeEnReIndicVec and rGeEnReIndicVec that
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
if nargin < 2 || isempty(enzymes)
    if ~isfield(cnap,'enzymes')
        [cnap,enzymes] = CNAgenerateGERassociation(cnap);
    else
        enzymes = cnap.enzymes;
    end
end
numEnzymes = length(enzymes);
if nargin < 3 || isempty(flagGeneOrEnzymeKO) || flagGeneOrEnzymeKO~='e'
    flagGeneOrEnzymeKO = 'g';
end
if nargin < 4
    koCost = ones(cnap.numr,1);
end
if nargin < 5
    kiCost = ones(cnap.numr,1);
end
if nargin < 6
    D = {};
    T = {};
end
if nargin < 8 || isempty(gkoCost) % number of genes is defined by 
    numGenes = max([enzymes(:).genes]);
    gkoCost = ones(numGenes,1);
    gkoCost_provided = 0;
else
    numGenes = length(gkoCost);
    gkoCost_provided = 1;  
end
if nargin < 9 || isempty(gkiCost)
    gkiCost = nan(numGenes,1);
    gkiCost_provided = 0;   
else, gkiCost_provided = 1;
end
if nargin < 10
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

if ~expand_notknockable
    affected_reacs = [];
    if ~isempty(gkoCost) % remove enzymes with notknockable genes
        for i = find(isnan(gkoCost(:)'))
            for j = 1:length(enzymes)
                if ismember(i,enzymes(j).genes)
                    affected_reacs = [affected_reacs enzymes(j).reactions];
                end
            end
        end
        for i = unique(affected_reacs)
            if sum(affected_reacs == i) == sum([enzymes(:).reactions] == i)
                koCost(i) = nan;
                kiCost(i) = nan;
            end
        end
        enzymes_w_noKO_genes = cellfun(@(x) all(ismember(find(isnan(gkoCost)),x)), {enzymes(:).genes});
    else
        enzymes_w_noKO_genes = zeros(1,length(enzymes));
    end
    % remove enzymes with only notknockable reactions
    enzymes_w_noKO_reacs = cellfun(@(x) all(ismember(find(isnan(koCost)),x)), {enzymes(:).reactions});
    enzymes = enzymes(  ~(enzymes_w_noKO_genes  | enzymes_w_noKO_reacs)  );
end

if isempty(enzymes)
    error('Model extension with enzymes and genes not possible. Make sure that at least one enzyme can be knocked out.');
end

% set gene and enzyme names
if ~isfield(enzymes,'strGene')
    geneNames  = strcat('G',num2str(1:numGenes));
else
    [a,b] = unique([enzymes(:).genes]);
    gname = [enzymes(:).strGene]';
    geneNames = arrayfun(@(x) ['gene ' num2str(x)],1:numGenes,'UniformOutput',0);
    geneNames(a) = gname(b);
end
if ~isfield(enzymes,'name')
    enzNames = strcat('E',num2str(1:length(enzymes)));
else
    enzNames = [enzymes(:).name]';
    if any(ismember(enzNames,geneNames))
        enzNames = strcat('Enz-',enzNames);
    end
end

% duplicate reactions that have an enzyme association and are reversible or 
% with isoenzyme association.
 % number of enzymes that catalyze a reaction:
enz_per_reac = hist([enzymes(:).reactions],1:cnap.numr);
 % subset of these reactions that are reversible:
rev_reac_w_enz = enz_per_reac & (sign(cnap.reacMin) .*  sign(cnap.reacMax) == -1)';
 % which reactions need to be split, create mapping. Which reactions are enzymatically catalyzed:
r_old2new = num2cell(1:cnap.numr);
r_old_enz = num2cell(zeros(1,cnap.numr));
for i = find(enz_per_reac) % indicate which reactions are cloned (+/-)
    r_old2new{i} = [i*ones(1,enz_per_reac(i)) -i*ones(1,enz_per_reac(i)*rev_reac_w_enz(i))]; % reactions
    corresp_enz_idx = find(arrayfun(@(x) ismember(i,x.reactions),enzymes));
    r_old_enz{i} = repmat(corresp_enz_idx,1,rev_reac_w_enz(i)+1); % new enzyme mapping
end
r_new2old = [r_old2new{:}];
r_new_enz = [r_old_enz{:}];
% reaction Type
rType(1:length(r_new2old)) = 'r'; % enzyme catalyzed reactions
% species Type
sType(1:cnap.nums) = 'm'; % metabolites

reacMap_old_new = full(sparse(abs(r_new2old),1:length(r_new2old),sign(r_new2old)));
 % reactions that stem back from reversible reactions with enzymes (and were split)
 % 0: reaction was not split, 1: positive sense reaction part, -1: reverse sens reaction part
r_was_rev = rev_reac_w_enz*reacMap_old_new;

if flagGeneOrEnzymeKO =='g' % GEM-Model
    gene_enz = num2cell(zeros(1,numGenes));
    for i = 1:numGenes
        gene_enz(i) = {find(cellfun(@(x) ismember(i,x),{enzymes(:).genes}))};
    end
    gene_reac = num2cell(zeros(1,numGenes));
    for i = 1:numGenes
        gene_reac(i) = {find(arrayfun(@(x) any(ismember(x,gene_enz{i})),r_new_enz))};
    end
end

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
gcnap.stoichMat = cnap.stoichMat* reacMap_old_new; % new Stoichmat
gcnap.reacMin = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMax = zeros(size(gcnap.stoichMat,2),1);
gcnap.reacMin(~r_was_rev) = cnap.reacMin(abs(r_new2old(~r_was_rev)));
gcnap.reacMin(~~r_was_rev) = 0;
gcnap.reacMax(r_was_rev >= 0) =  cnap.reacMax(abs(r_new2old(r_was_rev >= 0)));
gcnap.reacMax(r_was_rev <  0) = -cnap.reacMin(abs(r_new2old(r_was_rev <  0)));
gcnap.objFunc   = reacMap_old_new'*cnap.objFunc;
% reacIDs
% Counter for split reactions (only for naming)
countOcc = sum(cell2mat(arrayfun(@(x) cumsum(abs(r_new2old)==x).*(abs(r_new2old)==x),1:20,'UniformOutput',0)'),1)';
gcnap.reacNotes = cnap.reacNotes(abs(r_new2old));
gcnap.reacID    = cellstr(cnap.reacID(abs(r_new2old),:));
gcnap.reacID(r_new_enz>0,:) = strcat(gcnap.reacID(r_new_enz>0,:),'_enz_',num2str(countOcc(r_new_enz>0)));
gcnap.reacID(r_was_rev==-1,:) = strcat(gcnap.reacID(r_was_rev==-1,:),'_rev');
gcnap.reacID    = char(gcnap.reacID);
% specIDs
gcnap.specID    = cnap.specID;
gcnap.specNotes = cnap.specNotes;

%% 2. add enzymes and genes
%   extend stoichMat with m1 rows for enzymes (all have values -1 below N
%                                           and 1 to the bottom right of N)
%   extend stoichMat with m2 rows for genes (all bottom right corner)
%   extend stoichMat with n1 = m1 columns for gene-to-enzyme reactions
%   extend stoichMat with n2 = m2 columns for gene-producing reactions
%   entries in stoichMat as equivalent to enzyme 'consumption' through
%                                           reactions
% create vector that indicates whether reaction is gene expression, protein
% synthesis or reaction

% add enzmes
sType(end+1:end+numEnzymes) = 'e';
gcnap.stoichMat(sType=='e',:) = 0;
gcnap.specID = char([cellstr(gcnap.specID) ; enzNames(:)]);

% Add genes (if requested)
if flagGeneOrEnzymeKO =='g'
    sType(end+1:end+numGenes) = 'g';
    gcnap.stoichMat(sType=='g',:) = 0;
    gcnap.specID = char([cellstr(gcnap.specID) ; geneNames(:)]);
end

% add enzyme consumption to reactions
gcnap.stoichMat(sType=='e',:) = -sparse(r_new_enz(r_new_enz~=0),find(r_new_enz),1,numEnzymes,length(r_new_enz));

% add enzyme source terms / gene-enzyme terms
%  (fluxes entering the system are defined negative)
if flagGeneOrEnzymeKO =='g' % in GEM-Model
    rType(end+1:end+numEnzymes) = 'e';
    rType(end+1:end+numGenes)   = 'g';
    gcnap.stoichMat(sType=='e',rType=='e') = eye(numEnzymes); % add enzyme source pseudo-reaction
    gcnap.reacMin(rType=='e') = 0;
    gcnap.reacMax(rType=='e') = 1000;
    gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('GE-',enzNames)]);
    % change enzyme source terms into gene-enzyme pseudo-reactions
    gcnap.stoichMat(sType=='g',rType=='e') = -sparse(repelem(1:numGenes,cellfun(@length,gene_enz)),[gene_enz{:}],1,numGenes,numEnzymes);
    gcnap.stoichMat(sType=='g',rType=='g') = -eye(numGenes); % add gene source pseudo-reaction
    gcnap.reacMin(rType=='g') = -1000;
    gcnap.reacMax(rType=='g') = 0;
    gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('SG-',geneNames(:))]);
else % in EM-Model
    rType(end+1:end+numEnzymes) = 'e';
    gcnap.stoichMat(sType=='e',rType=='e') = -eye(numEnzymes); % add enzyme source pseudo-reaction
    gcnap.reacMin(rType=='e') = -1000;
    gcnap.reacMax(rType=='e') = 0;
    gcnap.reacID = char([cellstr(gcnap.reacID) ; strcat('SE-',enzNames)]);
end

gcnap.specNotes(size(gcnap.stoichMat,1)) = {''};
gcnap.reacNotes(size(gcnap.stoichMat,2)) = {''};
gcnap.objFunc(rType=='e' | rType=='g')   = 0;
reacMap_old_new = [reacMap_old_new, zeros(cnap.numr,numGenes+numEnzymes)];

%% 3. adapt kiCost, koCost, D, d, T, t
% kiCost and koCost vectors will be derived from reaction-Cost vectors
% if reactions have associated enzymes/genes, those will be marked as knockable.
% Else, the reactions will stay knockable. Conflicts are treated as follows:
% When geneA catalyzes R1 (knockable) and R2 (notknockable), the corresponding geneA 
% will be knockable. If R1 has a lower cost than R2, lower cost will be taken into account.

koCost_new = koCost(abs(r_new2old)); % translate koCost/kiCost to new set of reactions
kiCost_new = kiCost(abs(r_new2old));
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

for i = 1:size(gcnap.stoichMat,2)
    switch rType(i)
        case 'r' % reactions stay knockable if they were knockable before and don't have a gene-association
            if r_new_enz(i)==0
                gkoCost(i) = koCost_new(i);
                gkiCost(i) = kiCost_new(i);
            end
        case 'e'
            if flagGeneOrEnzymeKO =='e' % EM-Model
                enz_id = find(rType=='e')==i;
                enz_reac = find(ismember(abs(r_new2old),enzymes(enz_id).reactions));
                r_ko_costs = koCost_new(enz_reac);
                r_ki_costs = kiCost_new(enz_reac);
                if any(~isnan(r_ko_costs))
                    gkoCost(i) = min(r_ko_costs);
                end
                if any(~isnan(r_ko_costs))
                    gkiCost(i) = min(r_ki_costs);
                end
            end
        case 'g'
            gene_id = find(find(rType=='g') == i);
            if flagGeneOrEnzymeKO =='g' && ~ismember(gene_id,find(isnan(gkoCost))) % GEM-Model + gene is not already marked as notknockable
                r_ko_costs = koCost_new(gene_reac{gene_id});
                if any(~isnan(r_ko_costs)) && ~gkoCost_provided
                    gkoCost(i) = min(r_ko_costs);
                end
                r_ki_costs = kiCost_new(gene_reac{gene_id});
                if any(~isnan(r_ki_costs)) && ~gkiCost_provided
                    gkiCost(i) = min(r_ki_costs);
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
knockable_genes   = gkoCost(rType=='g');
knockable_enzymes = gkoCost(rType=='e');
not_knockable_warning = '';
for i = 1:cnap.numr % check if notknockable reactions still depend on at least one gene that is notknockable
    if isnan(koCost(i)) && any(r_old_enz{i})
        if flagGeneOrEnzymeKO =='e' % EM-Model
            if all(~isnan(knockable_enzymes(cellfun(@(x) any(ismember(i,x)),{enzymes(:).reactions}))))
                not_knockable_warning = [not_knockable_warning strtrim(cnap.reacID(i,:)) ...
                    ' :ENZ: ' num2str(r_old_enz{i}) ...
                    ' :GENE: ' num2str([enzymes(r_old_enz{i}).genes])  '; ' newline];
            end
        elseif flagGeneOrEnzymeKO =='g' % GEM-Model
            if all(~isnan(knockable_genes(cellfun(@(x) any(ismember(find(r_new2old==i),x)),gene_reac))))
                not_knockable_warning = [not_knockable_warning strtrim(cnap.reacID(i,:)) ...
                    ' :ENZ: ' num2str(r_old_enz{i}) ...
                    ' :GENE: ' num2str([enzymes(r_old_enz{i}).genes])  '; ' newline];
            end
        end
    end
end
if ~isempty(not_knockable_warning)
    warning(['A reaction that was before not knockable has become knockable through the expansion. '...
        'The affected reactions, enzymes and genes are ' newline not_knockable_warning 'Please check for ' ...
        'genetic overlaps of knockable and notknockable reactions before expanding the model.']);
end
% T and D
if ~isempty(T)
    for i = 1:length(T) % fill up to match reaction vector size
        gT{i}  = T{i}*reacMap_old_new;
    end
else
    gT = {};
end
if ~isempty(D)
    for i = 1:length(D) % fill up to match reaction vector size
        gD{i}  = D{i}*reacMap_old_new;
    end
else
    gD = {};
end

%% 4. in case of gui, generate map

if cnap.has_gui
    % deleting old boxes
    delete(unique(cnap.reacBoxes(:,4)))
    rBoxes = find(cnap.reacBoxes(:,5)>=1)';
    disp('redrawing reaction boxes');
    for reacIndex = rBoxes
        % declaration of handle
        cnan.open_projects.(cnap.net_var_name).gui.handles= struct;
        cnap = update_after_change(cnap);
        currmap = cnap.reacBoxes(reacIndex,5);
        % is visible/editable
        zw=cnap.reacBoxes(reacIndex,6);
        fig = cnap.figs(currmap,:);
        % make box
        zw1=uicontrol('Style', 'edit','Parent',fig(1), 'String', '###',...
            'Units','normalized','HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',cnap.reacID(reacIndex,:));
        set(zw1, 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
            {'check_for_right_click', 'reaceditmask'}, {'reacenr', reacIndex}});
        % save handle
        cnap.reacBoxes(reacIndex,4)=zw1;
        % adjust "zoom"
        place_box(fig,zw1,...
            cnap.reacBoxes(reacIndex,2),...
            cnap.reacBoxes(reacIndex,3),...
            cnap.reacFontSize(currmap),...
            cnap.reacBoxWidth(currmap),...
            cnap.reacBoxHeight(currmap));
        % put default rate in textbox
        if isnan(cnap.reacDefault(reacIndex))
            set(zw1,'String','#');
        else
            set(zw1,'String',num2str(cnap.reacDefault(reacIndex)));
        end
        % is visible/editable
        if(zw==2)
            set(zw1,'Style', 'text');
        elseif(zw==3)  %%non-visible
            set(zw1,'Visible','off');
        end
    end
end

% create gcnap Network
gcnap.epsilon = cnap.epsilon;
gcnap = CNAgenerateMFNetwork(gcnap);
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
