function p = plot_mcs_relationships(reac_names,mcs_1,mcs_2,mcs_3,mcs_4,mcs_5)
% Plot different sets of MCS side by side for easy comparison 

[~,~,~,a] = compare_mcs_sets(mcs_2,mcs_1);
[~,~,~,b] = compare_mcs_sets(mcs_3,mcs_1);
[~,~,~,c] = compare_mcs_sets(mcs_3,mcs_2);
e = ~(c*a) .* b;

% check for completeness
if all([max(a,[],2)>=2; 1]) && all([max(b,[],2)>=2; 1])
    disp('Hierarchy of MCS verified.');
% check for correctness
elseif all([max(a,[],2)~=1; 1]) && all([max(b,[],2)~=1; 1]) && all([max(c,[],2)~=1; 1])
    disp('Hierarchy of MCS verified, but some MCS must not have been found.');
else
    disp('No hierarchy of MCS.');
end

link_1_2            = cell(size(mcs_1,2),1);
link_1_2_type       = zeros(size(mcs_1,2),1);
[ro,co] = find(a);
for i = unique(co(:))'
    link_1_2{i} = ro(co==i)';
    link_1_2_type(i) = max(a(:,i));
end
link_2_3       = cell(size(mcs_2,2),1);
link_2_3_type  = zeros(size(mcs_2,2),1);
[ro,co] = find(c);
for i = unique(co(:))'
    link_2_3{i} = ro(co==i)';
    link_2_3_type(i) = max(c(:,i));
end
link_1_3         = cell(size(mcs_1,2),1);
link_1_3_type    = zeros(size(mcs_1,2),1);
[ro,co] = find(e);
for i = unique(co(:))'
    link_1_3{i} = ro(co==i)';
    link_1_3_type(i) = max(e(:,i));
end

% draw graph
G = digraph();
x = [];
y = [];
mcs_text = cell(size(mcs_1,2)+size(mcs_2,2)+size(mcs_3,2),1);
for i = 1:size(mcs_1,2)
    G = G.addnode(['MCS_set_1' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_weakGC(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_weakGC(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_1(:,i)))) ' (' mcs2text(reac_names,mcs_1(:,i)) ')'];
    x = [x 0];
    y = [y size(mcs_1,2)-i];
end
for i = 1:size(mcs_2,2)
    G = G.addnode(['MCS_set_2' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_directionalGC(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_directionalGC(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_2(:,i)))) ' (' mcs2text(reac_names,mcs_2(:,i)) ')'];
%     x = [x 3+i];
%     y = [y size(mcs_directionalGC,2)-i];
    x = [x 5];
    y = [y max(size(mcs_1,2),size(mcs_3,2))+size(mcs_1,2)-i];
end
for i = 1:size(mcs_3,2)
    G = G.addnode(['MCS_set_3' num2str(i)]);
%     mcs_text{G.numnodes,1} = num2str(sum(abs(mcs_substUp(:,i))));
%     mcs_text{G.numnodes,1} = mcs2text(reac_names,mcs_substUp(:,i));
    mcs_text{G.numnodes,1} = [num2str(sum(abs(mcs_3(:,i)))) ' (' mcs2text(reac_names,mcs_3(:,i)) ')'];
%     x = [x 3+i];
%     y = [y -3];
    x = [x 10];
    y = [y size(mcs_1,2)-i];
end
mcs_text = strrep(mcs_text,'EX\_o2\_e','O_2');

edges = nan(0,3); % source, target, type (identical = 3, subset = 2, superset = 1)
for i = 1:size(mcs_1,2)
    target = link_1_2{i};
    for j = target
        edges = [edges; i j+size(mcs_1,2) link_1_2_type(i)];
    end
    target = link_1_3{i};
    for j = target
        edges = [edges; i j+size(mcs_1,2)+size(mcs_2,2) link_1_3_type(i)];
    end
end
for i = 1:size(mcs_2,2)
    target = link_2_3{i};
    for j = target
        edges = [edges; i+size(mcs_1,2) j+size(mcs_1,2)+size(mcs_2,2) link_2_3_type(i)];
    end
end
color = nan(size(edges,1),1);
color(edges(:,3) == 1) = 8;
color(edges(:,3) == 2) = 1.5;
color(edges(:,3) == 3) = 0;
for i = 1:size(edges,1)
    if edges(:,3) ~= 1
        G = G.addedge(edges(i,1),edges(i,2));
    else
        G = G.addedge(edges(i,2),edges(i,1));
    end
end

% plot
p = plot(G,'XData',x,'YData',y,'EdgeCData',color,'NodeLabel',mcs_text);

% title
plot_title = {};
if ~isempty(mcs_1)
    plot_title = [plot_title {'MCS set 1'}];
end
if ~isempty(mcs_2)
    plot_title = [plot_title {'MCS set 2'}];
end
if ~isempty(mcs_3)
    plot_title = [plot_title {'MCS set 3'}];
end
title(strjoin(plot_title,'     vs     '));

% colors
% colormap('Turbo');
caxis([0 10]);

% layout
if ~isempty(mcs_1)
    sources = 1:size(mcs_1,2);
elseif ~isempty(mcs_2) && ~isempty(mcs_3)
    sources = 1:size(mcs_2,2);
else
    sources = [];
end
if ~isempty(mcs_3)
    sinks = size(mcs_1,2)+size(mcs_2,2)+(1:size(mcs_3,2));
elseif ~isempty(mcs_2) && ~isempty(mcs_1)
    sinks = size(mcs_1,2)+(1:size(mcs_2,2));
else
    sinks = [];
end
if ~isempty(sinks) && ~isempty(sources)
    layout(p,'layered','Direction','right','Sources',sources,'Sinks',sinks);
elseif isempty(sinks) && ~isempty(sources)
    layout(p,'layered','Direction','right','Sources',sources);
elseif ~isempty(sinks) && isempty(sources)
    layout(p,'layered','Direction','right','Sinks',sinks);
else
    layout(p,'layered','Direction','right');
end
% Prepare output and colors
p.NodeColor = 'black';
p.EdgeColor = 'black';
set(gca,'visible','off');

end

function ko_ki_text = mcs2text(reac_names,mcs)
    reac_names = strrep(cellstr(reac_names),'_','\_');
    kos = find(mcs<0);
    kis = find(mcs>0);
    if ~isempty(kis) 
        kis = [strjoin(strcat('+',reac_names(kis)),', ') ', '];
    else
        kis = [];
    end
    if ~isempty(kos) 
        kos = strjoin(reac_names(kos,:),', ');
    else
        kos = [];
    end
    ko_ki_text = [kis kos];
end
