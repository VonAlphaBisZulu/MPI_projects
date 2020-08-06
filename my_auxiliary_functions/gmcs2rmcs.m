function rmcs = gmcs2rmcs(gmcs,enzymes,rmap,rType)
%
% function 'gmcs2rmcs'
% ---------------------------------------------
% --> translates gene cut sets to reaction cut sets
%
% Usage: rmcs = gmcs2rmcs(gmcs,enzymes,rmap,rType)
%
% Input: 
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
% Output:
%
%  rmcs: provided as: num_reacs x num_mcs
% 

%
% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2020 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de
%
%% translate structure enzymes-array (n-1-m) to gpr_rules-array (1-1-m) if necessary
if isfield(enzymes,'reactions')
    gpr_rules = repelem_loc(enzymes,cellfun(@length,{enzymes.reactions}));
    count = repelem_loc(1:length(enzymes),cellfun(@length,{enzymes.reactions}));
    occ_count = sum(cell2mat(arrayfun(@(x) cumsum(abs(count)==x).*(abs(count)==x),unique(count),'UniformOutput',0)'),1)';
    for i = 1:length(gpr_rules)
        gpr_rules(i).reaction = gpr_rules(i).reactions(occ_count(i));
    end
    gpr_rules = rmfield(gpr_rules,'reactions');
else
    gpr_rules = enzymes;
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
    rule_ki_ko = zeros(1,length(gpr_rules));
    for j = 1:length(gpr_rules)
        if any(gmcs(gpr_rules(j).genes,i) == -1) && ~any(gmcs(gpr_rules(j).genes,i) == 1) && all(~isnan(gmcs(gpr_rules(j).genes,i)))
            rule_ki_ko(j) = -1;
        elseif all(gmcs(gpr_rules(j).genes,i) == 1) && ~any(isnan(gmcs(gpr_rules(j).genes,i)))
            rule_ki_ko(j) = 1;
        elseif any(isnan(gmcs(gpr_rules(j).genes,i)))
            rule_ki_ko(j) = nan;
        end
    end
    % trace back knocked-out or knocked-in reactions
    rlist = [gpr_rules(:).reaction];
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

% local implementation of the MATLAB repelem function.
function U = repelem_loc(V,varargin)
% This function is used to guarantee compatibility with MATLAB2014b and below
    % if V = 1xn or nx1 and N = 1xn
    if size(V,1) == 1 && length(varargin) == 1
        varargin{2} = 1;
        varargin = flip(varargin);
    elseif size(V,2) == 1 && length(varargin) == 1
        varargin{2} = 1;
    end
    % if V = 1xn and N = 1x1
    for i = find(cellfun(@length,varargin) == 1)
        varargin{i} = varargin{i}*ones(1,getelements(size(V)',i));
    end
	[reps{1:numel(varargin)}] = ndgrid(varargin{:});
	U = cell(cellfun(@length,varargin));
	U(:) = arrayfun(@(x) V(x)*ones(cellfun(@(y) y(x),reps)) ,1:numel(V),'UniformOutput',false);
	U = cell2mat(U);
end