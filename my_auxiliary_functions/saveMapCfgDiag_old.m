function saveMapCfgDiag( cnap )
% Usage:  saveMapCfgDiag( cnap )
% Input:
%   cnap:     CNA project object
%
% Opens a dialog box that asks to save the current map configuration
% permanently. If yes is chosen, the function calls
% CNAsaveProjectWithNewMap(cnap).
%
% Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
% Jun 13, 2017

global cnan;
try
    close(cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag);
    cnan.open_projects = rmfield(cnan.open_projects,[cnap.net_var_name,'_dialogs']);
catch
end
cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag = ...
    figure('Position',[800 500 460 150],'Name','Save Map configuration permanently',...
    'NumberTitle','off','Menubar','none','Resize','off');
        
        txt = uicontrol('Parent',cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag,...
            'Style','text',...
            'Position',[30 60 420 60],...
            'String','Maps were edited, added or removed from the project. Do you want to save map changes permanently? If yes, please also save project changes later.');
        
        btn1 = uicontrol('Parent',cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag,...
            'Position',[20 20 120 30],...
            'String','Save',...
            'Callback',strcat(cnap.net_var_name,' = CNAsaveProjectWithNewMap(',cnap.net_var_name,');'));
        
        btn2 = uicontrol('Parent',cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag,...
            'Position',[140 20 180 30],...
            'String','Save & Close Dialog',...
            'Callback',[cnap.net_var_name,' = CNAsaveProjectWithNewMap(',cnap.net_var_name,');',...
            'close(cnan.open_projects.',cnap.net_var_name,'_dialogs.saveMapDiag);',...
            'cnan.open_projects = rmfield(cnan.open_projects,''',cnap.net_var_name,'_dialogs',''');']);

        btn3 = uicontrol('Parent',cnan.open_projects.([cnap.net_var_name,'_dialogs']).saveMapDiag,...
            'Position',[320 20 120 30],...
            'String','Cancel',...
            'Callback',['close(cnan.open_projects.',cnap.net_var_name,'_dialogs.saveMapDiag);',...
            'cnan.open_projects = rmfield(cnan.open_projects,''',cnap.net_var_name,'_dialogs',''');']);
end

