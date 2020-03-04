cnap = CNAloadNetwork({'multiSsmall';1},1,1);

p_ex = findStrPos(cnap.reacID,'rp_ex');
t_ex = findStrPos(cnap.reacID,'rt_ex');
t_up = findStrPos(cnap.reacID,'rt_up');
s_up = findStrPos(cnap.reacID,'rs_up');
bm   = findStrPos(cnap.reacID,'rbm');
z    = findStrPos(cnap.reacID,'rz');

ft = 1;
fs = 1;
Y  = 1;

r{1} = findStrPos(cnap.reacID,'r[3,6]','regex');
r{2} = findStrPos(cnap.reacID,'r[3,6]|t_up','regex');
r{3} = findStrPos(cnap.reacID,'r[1,3]','regex');
r{4} = findStrPos(cnap.reacID,'r[1,3]|t_ex','regex');
% cnap.reacMin(r) = 0;
% cnap.reacMax(r) = 0;
for i = 1:4
%     cnapn = cnap;
%     cnapn.reacMin(r{i}) = 0.001;
%     cnapn.reacMax(r{i}) = 0.001;
    fixedFluxes = nan(1,cnap.numr);
    fixedFluxes(r{i}) = 0;

    c = full(sparse( 1 , p_ex, -1 ,1,cnap.numr));
    d = full(sparse( 1 , [t_up, s_up], [Y/ft Y/fs] ,1,cnap.numr));

    [maxyield,flux_vec] = CNAoptimizeYield(cnap,c,d,fixedFluxes,[],2,-1);
    disp(-maxyield);
end