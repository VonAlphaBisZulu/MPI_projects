function cnap = CNAloadMCSinGUI(cnap,mcs)
%% sets the entry in all reaction boxes for reactions that are defined in 
%  an mcs vector to 0 or a custom value (such as #)
    cnap.local.cutsets=mcs';
    cnap.local.num_cutsets=size(mcs,2);
    assignin('base',cnap.net_var_name,cnap);
    evalin('base',[cnap.net_var_name ' = show_cutsets(' cnap.net_var_name ');']);
end