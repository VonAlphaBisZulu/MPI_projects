weight = ...
[0,	741,	630,	406,	325,	325,	351,	0,	300,	276;
741,	0,	406,	351,	276,	0,	0,	0,	0,	0;
630,	406,	0,	300,	0,	276,	300,	0,	0,	0;
406,	351,	300,	0,	0,	0,	0,	0,	0,	0;
325,	276,	0,	0,	0,	0,	0,	0,	0,	0;
325,	0,	276,	0,	0,	0,	276,	0,	0,	0;
351,	0,	300,	0,	0,	276,	0,	0,	0,	0;
0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
300,	0,	0,	0,	0,	0,	0,	0,	0,	0;
276,	0,	0,	0,	0,	0,	0,	0,	0,	0];
counts = [56,42,50,34,32,31,30,29,28,28];
nnodes = length(counts);

screen_names = cellfun(@num2str,num2cell(1:10),'UniformOutput',false);
G = digraph(weight);
LWidths = 10*G.Edges.Weight/max(G.Edges.Weight);
h = plot(G,'NodeLabel',cellstr(screen_names),'LineWidth',LWidths);
for i = 1:nnodes
    highlight(h,i, 'MarkerSize', counts(i))
end