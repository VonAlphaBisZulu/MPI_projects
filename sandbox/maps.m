model_name = 'ivnat'; % iab or ivnat
%% remove old stuff
openf = fopen('all');
arrayfun(@fclose,openf);
fileID = fopen('networks','r');
text = char(fread(fileID)');
fclose(fileID);
text = regexprep(text,['^' model_name '.*\n'],'');
fileID = fopen('networks','w');
fprintf(fileID,text);
fclose(fileID);
try
    rmdir(model_name,'s');
end
clearvars -except model_name
clc;
startcna(1);

%% 3 cases:
%% 1. extend existing model with reactions from xls-file
CNAloadNetwork({'iJO1366';1});
waitfor(msgbox('Add reactions from .xls file'));
iJO1366 = CNAaddSpecsAndReacsFromFile(iJO1366,which('P14_Isobutanol.xls'));
waitfor(msgbox('remove reactions'));
iJO1366 = CNAremoveSpecsAndReacsDefinedInFile(iJO1366,which('P14_Isobutanol.xls'));
waitfor(msgbox('close project (don''t save)'));
closeDialog2('UserData','iJO1366');

%% 2. remap reactions to newly generated map (also the case when a large number of reactions is added to a network)
CNAloadNetwork({'iJO1366';1});
waitfor(msgbox('remap reactions to newly generated map'));
iJO1366.reacBoxes(300:700,5) = -1;
CNAgenerateMap(iJO1366);
waitfor(msgbox('close project (don''t save)'));
closeDialog2('UserData','iJO1366');

%% 3. build model from sbml
waitfor(msgbox('load sbml network and generate maps. Then save network'));
cnap = CNAsbmlModel2MFNetwork(which('iVnat.xml')); % load SBML iVnat or iAB
cnap.path = model_name;
CNAsaveNetwork(cnap); % Save model without gui
cnap = CNAgenerateMap(cnap); % Generate reaction map
closeDialog2('UserData',model_name);
waitfor(msgbox('add new project to networks'));
fileID = fopen('networks','a');
fprintf(fileID,[newline model_name char(9) model_name char(9) '1']);
fclose(fileID);
startcna;