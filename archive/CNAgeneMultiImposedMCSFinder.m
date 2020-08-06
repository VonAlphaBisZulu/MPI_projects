function [full_mcs, full_cnap, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, status, obj] = ...
    CNAgeneMultiImposedMCSFinder(cnap,T,t,D,d,koCost,kiCost,maxNumInterv,timeLimit,maxSolutions,use_compression,use_bigM,enum_method,gkoCost,gkiCost,gr_rules,verbose)
if nargin < 11 || isempty(use_compression)
    use_compression = 2;
end
if nargin < 12 || isempty(use_bigM)
    use_bigM = 0;
end
if nargin < 13 || isempty(enum_method)
    enum_method = 1;
end
if nargin < 17 || isempty(verbose)
    verbose = 1;
end
if isempty(kiCost)
    kiCost = nan(cnap.numr,1);
end
if isempty(D)
    lb_D={};
    ub_D={};
end
disp('== gene MCS Computation ==');
if nargin < 16 || isempty(gr_rules)
    if ~isfield(cnap,'gr_rules') || ~isfield(cnap,'genes')
        [~,~,genes,gr_rules] = CNAgenerateGERassociation(cnap);
    else
        gr_rules = cnap.gr_rules;
        genes = cnap.genes;
    end
else
    [a,b] = unique([gr_rules(:).genes]);
    gname = [gr_rules(:).strGene]';
    genes(a) = gname(b);
end
if nargin < 14 || isempty(gkoCost) % derive from reaction knock-outs, if not defined. Set reactions with gene-association to notknockable.
    gkoCost = ones(1,length(genes));
    for i = 1:length(genes) % if all associated reactions are notknockable, so is the gene.
        rules_1 = cellfun(@(x) ismember(i,x),{gr_rules(:).genes});  % bool vector: Rules that contain the gene
        if all(isnan(koCost([gr_rules(rules_1).reaction])))
            gkoCost(i) = nan;
        end
    end
else
    gkoCost = gkoCost(:)';
end
if nargin < 15 || isempty(gkiCost) % derive from reaction knock-ins, if not defined
    gkiCost = nan(1,length(genes));
    for i = 1:length(genes)
        rules_1 = cellfun(@(x) ismember(i,x),{gr_rules(:).genes});  % bool vector: Rules that contain the gene
        if ~all(isnan(kiCost([gr_rules(rules_1).reaction])))
            gkiCost(i) = min(kiCost([gr_rules(rules_1).reaction]));
        end
    end
else
    gkiCost = gkiCost(:)';
end
%% Prepare koCost, gkoCost, kiCost and gkiCost
gkoCost(~isnan(gkiCost)) = nan; % gene knock-ins 'override' gene knock-outs
c_macro = cnap.macroDefault;
% Transform ki & ko vector to row vector with 'nan' as notknockable and a double value as 'weight'
kiCost = kiCost(:)';
kiCost([gr_rules(:).reaction]) = nan;  % gene-associations 'override' reaction knock-outs
koCost = koCost(:)';
koCost([gr_rules(:).reaction]) = nan;
koCost(~isnan(kiCost)) = nan; % knock-ins 'override' knock-outs

% Making backups of some variables. Needed later for mcs expansion
full_koCost = koCost;
full_kiCost = kiCost;
full_gkoCost = gkoCost;
full_gkiCost = gkiCost;
full_T = T;
full_D = D;
cnap_orig = cnap;
gr_rules_orig = gr_rules;

%% Determine boundaries of metabolic model
disp('FVA to determine model bounds.');
[cnap.reacMin, cnap.reacMax] = CNAfluxVariability(cnap,[],c_macro,2);

%% Determine boundaries in desired scenarios and identify further essentials/notknockables
disp('FVA to determine model bounds under desired constraints.');
essential = zeros(1,cnap.numr);
for i = 1:length(D)
    [lb_D{i}, ub_D{i}] = CNAfluxVariability(cnap,[],c_macro,2,1:cnap.numr,D{i},d{i}); %#ok<*AGROW> Allow change of size in each loop iteration
    lb_D{i}( abs(lb_D{i})<=cnap.epsilon ) = 0;
    ub_D{i}( abs(ub_D{i})<=cnap.epsilon ) = 0;
    essential = (sign(lb_D{i}).*sign(ub_D{i})) == 1;
    koCost(essential) = nan; % make essential reactions "notknockable"
    % *could also be used to select mandatory knock-ins
end

%% 1. Reduce and lump set of rules
disp('Reducing set of gene reaction rules.');
rule_removelist = zeros(1,length(gr_rules));
[reac_abund,reac] = hist([gr_rules(:).reaction],unique([gr_rules(:).reaction]));
cgenes = {gr_rules(:).genes}; % cell array of genes
for i = 1:length(genes)
    % if a gene only catalyzes essential reactions -> set gene notknockable
    rules_1 = cellfun(@(x) ismember(i,x),cgenes); % bool vector: Rules that contain the gene
    if all(essential([gr_rules(rules_1).reaction]))
        gkoCost(i) = nan;
    end
    % if the gene is essential to one essential reaction -> set gene notknockable
    for j = unique([gr_rules(rules_1).reaction])
                         % if all reaction-gene-rules for a particular reaction contain the gene
        if essential(j) && sum([gr_rules(rules_1).reaction] == j) == reac_abund(j == reac)
            gkoCost(i) = nan;
        end
    end
    % if gene is notknockable remove it from all rules
    if isnan(gkoCost(i)) && isnan(gkiCost(i)) && any(rules_1)
        for j = find(rules_1)
            gr_rules(j).genes = setdiff(gr_rules(j).genes,i);
        end
    end
end
% if all genes of one rule are notknockable, delete all rules for the same reaction
% because the reaction can never be knocked out.
% Here, through preprocessing, some rules don't contain any genes anymore (genes = []) 
% and are thus also notknockable
for i = 1:length(gr_rules)
    if all(isnan(gkoCost(gr_rules(i).genes))) && all(isnan(gkiCost(gr_rules(i).genes))) % finds also 'empty' (genes = []) rules 
        rule_removelist([gr_rules(:).reaction] == gr_rules(i).reaction) = 1;
    end
    % remove rules for reactions that are off
    if cnap.reacMin(gr_rules(i).reaction) == 0 && cnap.reacMax(gr_rules(i).reaction) == 0
        rule_removelist(i) = 1;
    end
end
gr_rules = gr_rules(~rule_removelist);

% genes that don't occur in rules are notknockable
gkoCost(setdiff(1:length(gkoCost),[gr_rules(:).genes])) = nan;
gkiCost(setdiff(1:length(gkiCost),[gr_rules(:).genes])) = nan;

%% 2. Join genes, adapt reaction's ko- and kiCosts
if ~isempty(gr_rules)
    gr_rules = remove_nonminimal(gr_rules); 
    [gene_subst_mat_AND,gr_rules,gkoCost,gkiCost] = unite_genes_AND(gr_rules,gkoCost,gkiCost);
    [gene_subst_mat_AND_2,gene_subst_mat_OR, gr_rules,gkoCost,gkiCost] = unite_genes_OR( gr_rules,gkoCost,gkiCost,cnap);
else
    gene_subst_mat_AND = eye(length(genes));
    gene_subst_mat_AND_2 = eye(length(genes));
end

if use_compression > 1
%% compress network if indicated
    disp('Compressing GEM model.');
    [        cmp_cnap, T, D, koCost, kiCost, cmp_transf_mat, lb_D, ub_D, gkoCost, gkiCost, gr_rules ] = ... Compression
    compress(cnap, T, D, koCost, kiCost,                 lb_D, ub_D, gkoCost, gkiCost, gr_rules);
    text = 'compressed ';
else
    text = '';
    cmp_transf_mat = eye(cnap.numr);
    cmp_cnap = cnap;
end

%% test again feasibility of target and desired vectors
disp(['Verifying D and T region feasibility in ' text 'GEM model.']);
testRegionFeas(cmp_cnap,c_macro,T,t,D,d);

%% incorporate genes and enzymes
disp('Generating GEM model.');
if ~isempty(gr_rules)
    [ cmp_cnap, M_GEM_map, ~,koCost, kiCost, T, D, rType ] = ...
    CNAextendModelGenes( cmp_cnap, gr_rules, koCost, kiCost, T, D, gkoCost, gkiCost, 0);
    for i = 1:length(D)
        [lb_D{i}, ub_D{i}] = tailorBounds(lb_D{i},ub_D{i},M_GEM_map);
    end
else
    M_GEM_map = eye(cmp_cnap.numr);
    rType = repmat('r',1,cmp_cnap.numr);
end

%% compute MCS
[cmp_mcs, status, obj] = CNAmultiImposedConstrISfinder(cmp_cnap,T,t,D,d,koCost,kiCost,maxNumInterv,timeLimit,maxSolutions,use_compression,use_bigM,enum_method,verbose);

if status ~= 2
    cmp_mcs_exp = sparse(cmp_mcs);
    cmp_mcs_exp(isnan(cmp_mcs)) = 0; % for expansion non-knock-ins are first set to 0

    cmp_cnap.mcs.kiCost = kiCost;
    cmp_cnap.mcs.koCost = koCost;
    cmp_cnap.mcs.T = T;
    cmp_cnap.mcs.t = t;
    cmp_cnap.mcs.D = D;
    cmp_cnap.mcs.d = d;

    disp('Expand MCS.');
    %% expand to gene rules
    mcs_idx_cmp_full = 1:size(cmp_mcs_exp,2);
    % 1st undo gene extension and compression (no expansion, because knockable reactions without genes were not compressed)
    mcs_reac  = sparse(cmp_transf_mat*M_GEM_map*cmp_mcs_exp);
    mcs_genes = sparse(cmp_mcs_exp(rType == 'g',:));
    % prepare expansion:
    % replacements
    % OR -> KI are expanded, AND -> KO are substituted
    %                             number of possible ways to employ functional knock-out/in
    % substitution array {    1    x      {  1     x     n  } }
    %                                             genes needed for one particular way
    gene_subst_mat_AND_combined = logical(double(gene_subst_mat_AND)*double(gene_subst_mat_AND_2));
    subst = cell(1,sum(rType == 'g'));
    for i = 1:sum(rType == 'g')
        if ~isnan(kiCost(sum(rType ~= 'g')+i)) % KI
            % find rule: (X OR Y OR Z), KI of only one needed for guaranteeing function
            subst{i} = num2cell(  find(gene_subst_mat_OR(:,i))  )'; % KI:OR
            for j = 1:size(subst{i},2)
                subst{i}{j} = find(gene_subst_mat_AND_combined(:,subst{i}{j}))'; % KI:AND
            end
        elseif ~isnan(koCost(sum(rType ~= 'g')+i)) % KO
            % find rule: (X OR Y OR Z), KO of all needed for blocking function
            subst(i) = { find(gene_subst_mat_OR(:,i))' }; % KO:OR
            % find all genes (Xa AND Xb AND Xc) for each sub-rule (X OR Y OR Z) and explore all possible
            % combinations to falsify the rule
            summands = cell(1,size(subst{i},2));
            for j = 1:size(subst{i},2)
                summands{j} = find(gene_subst_mat_AND_combined(:,subst{i}(j)))'; % KO:AND
            end
            combinations = combvec(summands{:})';
            combinations = arrayfun(@(x) unique(combinations(x,:)),1:size(combinations,1),'UniformOutput',false); % eliminate reduncancies
            isminimal = true(1,length(combinations));
            for j = 1:length(combinations) % delete non-minimals
                for k = 1:(j-1)
                    if all(ismember(combinations{k},combinations{j}))
                        isminimal(j) = false;
                    end
                end
            end
            subst{i} = combinations(isminimal);
        end
    end
    % expand mcs
    mcs_genes2 = sparse([],[],[],size(gene_subst_mat_AND_combined,1),size(mcs_genes,2));
    for i = find(any(mcs_genes,2))'
        num_split = length(subst{i})*(mcs_genes(i,:)~=0) + double(mcs_genes(i,:)==0); % how many times each mcs is copied
        a = repelem((1:size(mcs_genes,2)) .* double(mcs_genes(i,:)~=0),num_split); % map of reactions before and after splitting
        occ = occurcount( full(a) ) .* double(a'~=0); % counter for each mcs (if split 3 times: (R1)1, (R1)2, (R1)3)
        mcs_genes = repelem(mcs_genes,1,num_split); % make copies of mcs that need to be split
        mcs_genes2 = repelem(mcs_genes2,1,num_split); % make copies of mcs that need to be split
        for j = setdiff(unique(occ)',0)
            mcs_genes2(subst{i}{j}, occ == j) = repmat(mcs_genes(i, occ == j),size(subst{i}{j},2),1);
        end
        mcs_idx_cmp_full = repelem(mcs_idx_cmp_full,num_split);
    end
    mcs_reac2 = mcs_reac(:,mcs_idx_cmp_full); % copy also mcs of reaction part to fit with the gene part

    %% Preparing output and eliminating MCS that are too expensive
    disp('Generating full model and preparing output.');
    if ~isempty(gr_rules)
        [ full_cnap, rmap, gmap, ~,~, full_T, full_D, rType ] = CNAextendModelGenes( cnap_orig, gr_rules_orig, [], [], full_T, full_D );
    else
        full_cnap = cnap_orig;
        rmap = eye(full_cnap.numr);
        gmap = [];
        full_T = T;
        full_D = D;
        rType = repmat('r',1,full_cnap.numr);
    end
        
    mcs_expanded = [rmap(:,rType=='r')'*mcs_reac2; zeros(sum(rType=='p'),size(mcs_reac2,2)) ; mcs_genes2];

    % translate kiCost and koCost
    [a,b] = find(rmap);
    full_kiCost(b) = full_kiCost(a);
    full_koCost(b) = full_koCost(a);
    ivCost = [full_kiCost nan(1,sum(rType=='p')) full_gkiCost];

    ivCost(~isnan(full_koCost))  = full_koCost(~isnan(full_koCost));
    ivCost([false(1,sum(rType~='g')) ~isnan(full_gkoCost)]) = full_gkoCost(~isnan(full_gkoCost));

    if maxNumInterv == inf
        maxNumInterv = 1e9; % workaround
    end
    ivCost(isnan(ivCost))    = maxNumInterv+1; % higher number, so that mcs containing this intervention are sortet out
    % check if knock-out and knock-in costs are still met
    mcs_affordable = sum((mcs_expanded~=0 & ~isnan(mcs_expanded)).*repmat(ivCost',1,size(mcs_expanded,2)),1) <= maxNumInterv;
    for i = find(~isnan(full_kiCost)) % when KI-reactions were not knocked-in, set them to NaN in output
        mcs_expanded(i,mcs_expanded(i,:) == 0) = nan;
    end
    full_mcs = mcs_expanded(:,mcs_affordable);
    mcs_idx_cmp_full = mcs_idx_cmp_full(mcs_affordable);
    disp(['MCS found: ' num2str(size(full_mcs,2)) '. (' num2str(size(cmp_mcs_exp,2)) ' compressed)']);
    % prepare output
    full_cnap.mcs.kiCost = full_kiCost;
    full_cnap.mcs.koCost = full_koCost;
    full_cnap.mcs.ivCost = ivCost;
    full_cnap.mcs.rmap = rmap;
    full_cnap.mcs.gmap = gmap;
    full_cnap.mcs.T = full_T;
    full_cnap.mcs.t = t;
    full_cnap.mcs.D = full_D;
    full_cnap.mcs.d = d;
else
    disp('no MCS found');
    full_mcs = nan;
    full_cnap = nan;
    cmp_mcs = nan;
    cmp_cnap = nan;
    mcs_idx_cmp_full = nan;
end
end

%% 1. compress
function [cmp_cnap, cmp_T, cmp_D, cmp_koCost, cmp_kiCost, cmp_mapReac, lb_D, ub_D, gkoCost, gkiCost, gr_rules] = compress(cnap,T,D,koCost,kiCost, lb_D, ub_D, gkoCost, gkiCost, gr_rules)
    %% Prepare / Protect some reactions from compression
    % identify essential reactions and adapt notknockable vector
    non_compress_reacs = any([cell2mat(D') ; cell2mat(T')],1);
    non_compress_reacs(~isnan(kiCost)) = true; % don't compress knock-in-able
    % don't compress knockable genes without gene affiliation (e.g. O2 uptake)
    non_compress_reacs(setdiff(find(~isnan(koCost)),[gr_rules(:).reaction])) = true;
    non_compress_reacs = find(non_compress_reacs);

    r_off      = abs(cnap.reacMin)<=cnap.epsilon & ...
                 abs(cnap.reacMax)<=cnap.epsilon;

    javastderr= java.lang.System.err;
    java.lang.System.setErr(java.io.PrintStream('cplex_stderr.log'));
    %% Compress
    [~,~,cmp_mapReac,~,cmp_cnap] = CNAcompressMFNetwork(cnap,non_compress_reacs,[],1,0,1,r_off,1);
    java.lang.System.setErr(javastderr);
    % remap MCS-region vectors
    cmp_T = cellfun(@(x) x*cmp_mapReac,T,'UniformOutput',0);
    cmp_D = cellfun(@(x) x*cmp_mapReac,D,'UniformOutput',0);
    lumpedReacs = double(cmp_mapReac ~= 0);
    lumpedReacs(lumpedReacs == 0) = nan;
    cmp_kiCost = arrayfun(@(x) sum(kiCost(~isnan(lumpedReacs(:,x)))),1:size(lumpedReacs,2)); % for KI, all lumped reactions need to be knocked in
    cmp_koCost = min(lumpedReacs.*koCost'); % otherwise min value would be 0
    for i = 1:length(D)
        [lb_D{i}, ub_D{i}] = tailorBounds(lb_D{i},ub_D{i},cmp_mapReac);
    end
    %% Compress rules - get gene-reaction rules for old model
    % rule: cell array
    %      {  1      x     n{  1    x     m  }  }
    %                  n: number of terms connected by OR
    %                            m: indices of genes connected by AND 
    numgenes = length(gkoCost);
    rule = repmat({logical.empty(numgenes,0)},1,size(cmp_mapReac,1)); % empty rules. For each reaction
    for i = unique([gr_rules(:).reaction])
        for j = 1:length(gr_rules)
            if i == gr_rules(j).reaction
                rule{i} = [rule{i} full(sparse(gr_rules(j).genes,1,true,numgenes,1))]; % add rules
            end
        end
    end
    %% Compress enzymes - map to compressed model (rules of lumped reactions are joined with AND)
    emptyRules = cellfun(@isempty,rule)';
    ruleMat = repmat({[]},1,cmp_cnap.numr);
    for i = 1:cmp_cnap.numr
        reacs = find(logical(cmp_mapReac(:,i)) & ~emptyRules); % find reaction rules to join
        if ~isempty(reacs)
            ruleMat{i} = logical(rule{reacs(1)}); % start with reaction rule of first of the lumped reactions
            for j = 2:length(reacs) % multiply rule for every AND that occurrs (equivalent to expansion of (a+b)*c = a*c + b*c)
                ruleMat{i} = repelem(ruleMat{i},1,size(rule{reacs(j)},2)) | repmat(rule{reacs(j)},1,size(ruleMat{i},2));
                if size(rule{reacs(j)},2) > 1 || j == length(reacs) % reduce number of rules if possible, everytime their number increases through expansion
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
    [a,b] = unique([gr_rules(:).genes]);
    gname = [gr_rules(:).strGene]';
    gNames = repmat({''},1,length(gkoCost));
    gNames(a) = gname(b);
    gr_rules = gr_rules([]);
    c = 1;
    for i = find(~cellfun(@isempty,ruleMat))
        for j = 1:size(ruleMat{i},2)
            gr_rules(c).reaction = i;
            gr_rules(c).strReac = cellstr(cmp_cnap.reacID(i,:));
            gr_rules(c).genes = find(ruleMat{i}(:,j))';
            gr_rules(c).strGene = gNames(gr_rules(c).genes);
            gr_rules(c).name = {['Rule-' strjoin(gr_rules(c).strGene,'-') '-r' strtrim(cmp_cnap.reacID(i,:))]};
            c = c+1;
        end
    end
    % sort rules for genes (actually not necessay, but more concurrent with original order)
    [~,b] = sort(cellfun(@min,{gr_rules(:).genes}));
    gr_rules = gr_rules(b);
end
%% 2. generate reduced bounds
function [cmp_lb, cmp_ub] = tailorBounds(lb,ub,map)
    [cmp_lb, cmp_ub] = deal(nan(size(map,2),1));
    for k=1:size(map,2)
        forw=find(map(:,k)>0);
        revs=find(map(:,k)<0);
        if(~isempty(forw))
            cmp_lb(k) = max(lb(forw)./map(forw,k));
            cmp_ub(k) = min(ub(forw)./map(forw,k));
        end
        if(~isempty(revs))
            cmp_lb(k) = max([cmp_lb(k);ub(revs)./map(revs,k)]);
            cmp_ub(k) = min([cmp_ub(k);lb(revs)./map(revs,k)]);
        end
    end
end
%% 3.1 remove non-minimal rules
% A or (A and B) = A
function gr_rules = remove_nonminimal(gr_rules)
    maxreac = max([gr_rules(:).reaction]);
    reac_rule = sparse([gr_rules(:).reaction],1:length(gr_rules),1,maxreac,length(gr_rules));
    [~,~,rule_reac_abund] = unique(reac_rule','rows','stable'); 
    [rule_abund,isorule_map] = hist(rule_reac_abund,unique(rule_reac_abund));
    isorule_grule_ismin = false(1,length(gr_rules));
    % 1. for non-isoenzymes, the only existing gene rule is the shortest
    isorule_grule_ismin( ismember(rule_reac_abund', isorule_map(rule_abund == 1)') ) = 1;
    % 2. for isoenzymes, gene rules are compared to eliminate redundancies
    for i = isorule_map(rule_abund > 1)'
        genes_i = {gr_rules(rule_reac_abund == i).genes};
        isoenz = find(rule_reac_abund == i);
        isoenz_gene_rule = sparse(  cell2mat(genes_i),...
                                    repelem(1:length(isoenz), cellfun(@length,genes_i)),...
                                    1);
        % remove redundancies
        for j = 1:length(isoenz) % mark and eliminate redundand and non-minimal gene rules:
            if sum(all(repmat(isoenz_gene_rule(:,j),1,size(isoenz_gene_rule,2)) >= isoenz_gene_rule))==1  % A or (A and B) = A
                isorule_grule_ismin(isoenz(j)) = 1; % rule is minimal e.g. A
            else
                isorule_grule_ismin(isoenz(j)) = 0; % rule is not minimal e.g. (A and B)
            end
        end
    end
    gr_rules = gr_rules(isorule_grule_ismin); % remove redundant enzymes -> this also reduces set of genes
end
%% 3.2 unite and substitute genes if they occur always and only together with each other
% genes, connected with AND -> koCost are kept the same (lowest), kiCost are summed up
function [gene_subst_mat,gr_rules,gkoCost,gkiCost] = unite_genes_AND(gr_rules,gkoCost,gkiCost)
    genes_old = 1:length(gkoCost);
    [a,b] = unique([gr_rules(:).genes]);
    gname = [gr_rules(:).strGene]';
    gNames_old = repmat({''},1,length(genes_old));
    gNames_old(a) = gname(b);
    % remove redundant AND associations:
    gene_rule = cell2mat(cellfun(@(x) sparse(x,1,1,length(genes_old),1),   {gr_rules(:).genes}   ,'UniformOutput',0));
    gene_rule = logical(gene_rule);
    [a,~,gene_subst_mat] = unique(gene_rule,'rows','stable'); % unite equivalent genes
    gene_subst_mat = full(sparse(gene_subst_mat, 1:length(gene_subst_mat), true)); % new gene mapping    
    gene_subst_mat = gene_subst_mat(any(a,2),:); % keep genes with rule association
    genes_new = 1:size(gene_subst_mat,1); % genes_1: Genes-Vector after removing AND-redundancies
    gkoCost = arrayfun(@(x) min(gkoCost(gene_subst_mat(x,:))), genes_new); % take lowest cost for knock-out
    gkiCost = arrayfun(@(x) sum(gkiCost(gene_subst_mat(x,:))), genes_new); % take sum of costs for knock-in
    gNames_new =  arrayfun(@(x) ['(',strjoin(gNames_old(gene_subst_mat(x,:)),'_and_'),')'],genes_new,'UniformOutput',0);
    % reassign gene indices
    for i = 1:length(gr_rules)
        gr_rules(i).genes = unique(arrayfun(@(x) find(gene_subst_mat(:,x)), gr_rules(i).genes),'stable');
        gr_rules(i).strGene = gNames_new(gr_rules(i).genes);
    end
    gene_subst_mat = gene_subst_mat';
end
%% 3.3 unite and substitute genes if they occur only in one reaction
% genes, connected with OR -> koCost are summed up, kiCost are kept the same (lowest)
function [genes_to_join_AND,genes_to_join_OR,gr_rules,gkoCost,gkiCost] = unite_genes_OR(gr_rules,gkoCost,gkiCost,cnap)
    numGenes_old = length(gkoCost);
    numGenes     = numGenes_old;
    genes_old = 1:numGenes_old;
    [a,b] = unique([gr_rules(:).genes]);
    gname = [gr_rules(:).strGene]';
    gNames = repmat({''},1,length(genes_old));
    gNames(a) = gname(b);
    maxReac = max([gr_rules(:).reaction]);
    ruleMat = repmat({[]},1,maxReac);
    for i = unique([gr_rules(:).reaction])
        rules = [gr_rules(:).reaction] == i;
        rule_y = {gr_rules(rules).genes};
        rule_x = repelem(1:sum(rules),cellfun(@length,rule_y));
        ruleMat{i} = sparse([rule_y{:}],rule_x,true,numGenes_old,sum(rules));
    end
    ruleMat_all = [ruleMat{:}];
    ruleMat_new = ruleMat;
    gene_abund_all = sum(ruleMat_all,2);
    % extend genes and rules vector first, then remove obsolete entries
    genes_to_join_AND = logical(eye(numGenes_old));
    genes_to_join_OR  = logical(eye(numGenes_old));
    for i = unique([gr_rules(:).reaction]) % for all reactions with rules
        if size(ruleMat{i},2) > 1  % if reaction has multiple rules ('or')
            gene_abund_1 = sum(ruleMat{i},2);
            % check if all gene occurrencies are with the rules of this reaction. If yes, genes can
            % be joined.
            if gene_abund_1(logical(gene_abund_1)) == gene_abund_all(logical(gene_abund_1))
                for k = find(sum(ruleMat{i},1)>1) % first join genes that are connected with AND
                    % create new gene
                    numGenes = numGenes+1;
                    genes_to_join_AND( ruleMat{i}(:,k) ,numGenes) = 1;
                    gNames{numGenes} = ['[',strjoin(gNames(ruleMat{i}(:,k)),'_and_'),']'];
                    gkoCost(numGenes) = min(gkoCost(ruleMat{i}(:,k)));
                    gkiCost(numGenes) = sum(gkoCost(ruleMat{i}(:,k)));
                    % substitute genes in gene rule
                    ruleMat_new{i}(  :  ,k) = 0;
                    ruleMat_new{i}(numGenes,k) = 1;
                end
                % join genes that are connected with 'or'
                numGenes = numGenes+1;
                genes_to_join_OR( any(ruleMat_new{i},2), numGenes ) = 1;
                gNames{numGenes} = ['[',strjoin(gNames(any(ruleMat_new{i},2)),'_or_'),']'];
                ruleMat_new{i} = false(numGenes,1);
                ruleMat_new{i}(numGenes) = true;
                % calculate KO- or KI-costs
                originalRules  = genes_to_join_AND(:,genes_to_join_OR(:,ruleMat_new{i}));
                if all(~isnan(gkoCost(logical(gene_abund_1)))) % KO -> minimal hitting set problem with exhaustive search
                    gkiCost(numGenes) = nan;
                    occurringGenes = any(originalRules,2);
                    geneCombinations = false(2^sum(occurringGenes),numGenes);
                    geneCombinations(1:2^sum(occurringGenes),occurringGenes) = dec2bin(2^sum(occurringGenes)-1:-1:0)-'0';
                    cost = nan(size(geneCombinations,1),1);
                    for j = 1:size(geneCombinations,1)
                        if all(any(originalRules(geneCombinations(j,:),:),1)) % if all rules are hit
                            cost(j) = sum(gkoCost(geneCombinations(j,:)));
                        end
                    end
                    gkoCost(numGenes) = min(cost);
                elseif all(~isnan(gkiCost(logical(gene_abund_1)))) % KI -> minimal hitting set problem
                    gkoCost(numGenes) = nan;
                    cost = nan(size(originalRules,2),1);
                    for j = 1:size(originalRules,2)
                        cost(j) = sum(gkiCost(originalRules(:,j)));
                    end
                    gkiCost(numGenes) = min(cost);
                else
                    warning(['Ambigous knockability indication. Reaction ' num2str(i) ...
                             ' has genes that can be knocked in and others that can be knocked out.']);
                end
            end
        end
    end
    gr_rules = gr_rules([]);
    c = 1;
    for i = find(~cellfun(@isempty,ruleMat_new))
        for j = 1:size(ruleMat_new{i},2)
            gr_rules(c).reaction = i;
            gr_rules(c).strReac = cellstr(cnap.reacID(i,:));
            gr_rules(c).genes = find(ruleMat_new{i}(:,j))';
            gr_rules(c).strGene = gNames(gr_rules(c).genes);
            gr_rules(c).name = {['Rule-' strjoin(gr_rules(c).strGene,'-') '-r' strtrim(cnap.reacID(i,:))]};
            c = c+1;
        end
    end
    % sort rules for genes
    [~,b] = sort(cellfun(@min,{gr_rules(:).genes}));
    gr_rules = gr_rules(b);
    gkoCost(setdiff(1:length(gkoCost),[gr_rules(:).genes])) = nan; % genes that don't occur in enzymes anymore
    gkiCost(setdiff(1:length(gkiCost),[gr_rules(:).genes])) = nan; % are set to notknockable (and are not added to the model later).
end
%% 4. test region feasibility
function testRegionFeas(cnap,c_macro,T,t,D,d)
    for i = 1:length(t)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, T{i}, t{i}))
            error(['At least one target region (T' num2str(i) ') is infeasible in the original model']);
        end
    end
    for i = 1:length(d)
        if isnan(CNAoptimizeFlux(cnap, [], c_macro, 2, -1, 0, D{i}, d{i}))
            error(['At least one desired region (D' num2str(i) ') is infeasible in the original model']);
        end
    end
end
%% 5. write enzymes-vector back to model % not used, but might be helpful sometimes
function grRules = gen_grRules_from_enzymes(cnap,enzymes)
    grRules = cell(cnap.numr,1);
    [a,b] = unique([enzymes(:).genes]);
    gname = [enzymes(:).strGene]';
    gNames(a) = gname(b);
    for i = 1:cnap.numr
        rule = {};
        enz = find( cellfun(@(x) ismember(i,x),{enzymes.reactions}) );
        if ~isempty(enz)
            for j = enz
                str = strjoin(gNames(enzymes(j).genes),' and ');
                if length(enz) > 1 && length(enzymes(j).genes) > 1
                    str = ['(' str ')'];
                end
                rule = [rule str];
            end
            grRules{i} = strjoin(rule,' or ');
        else
            grRules{i} = '';
        end
    end
end
%% 6. resort gene indices after deleting enzymes
function [enzymes, gene_map] = rearrange_genes(enzymes)
    genes_old = unique([enzymes(:).genes]); % reassign gene indices
    genes_new(genes_old) = 1:length(genes_old);
    for i = 1:length(enzymes)
        enzymes(i).genes = genes_new(enzymes(i).genes);
    end
    gene_map = full(sparse(1:length(genes_old),genes_old,1,length(genes_old),length(genes_new)));
end
%% occurrency counter [1 2 2 3 3] -> [1 1 2 1 2]
function occ_count = occurcount(y)
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(y)==x).*(abs(y)==x),unique(y),'UniformOutput',0)'),1)';
end