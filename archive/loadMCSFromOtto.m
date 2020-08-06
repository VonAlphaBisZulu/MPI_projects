workingDir = pwd;
cd('/home/schneiderp/ottohome/Documents/MATLAB/CNA_SVN/_StrainBooster/_My_Simulations/Solutions')
files = dir;
mcsAndModelList = {};
for fidx = 1:length(files)
    [~,filename,ext] = fileparts(files(fidx).name);
    if strcmp(ext,'.mat')
        if ~strcmp(filename(end-4:end),'_cMCS')
            mcsAndModelList = [mcsAndModelList [filename ext]];
        end
    end
end
for fidx = 1:length(mcsAndModelList)
    loadedvars = load(char(mcsAndModelList(fidx)));
    tar = find(loadedvars.T(1,:)~=0);
    % get metabolite name
    stored(fidx).name = strrep(strrep(strtrim(loadedvars.iJO.reacID(tar(tar~=164),:)),'_LPAREN_e_RPAREN_',''),'R_EX_','');
    for fn = fieldnames(loadedvars)'
        stored(fidx).(char(fn)) = loadedvars.(char(fn));
    end
end
cd(workingDir)