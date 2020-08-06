% add Java to write xls
[jarpath,~,~] = fileparts(which('dom4j-1.6.1.jar'));
jf = dir(jarpath);
for i = 1:length(jf)
    if length(jf(i).name) >= 4
        if all(jf(i).name(end-3:end) == '.jar')
            javaaddpath([jarpath '/' jf(i).name]);
        end
    end
end

startcna(1)
cnap = CNAloadNetwork({'iJO1366';1},1,1);
cnapBU = cnap;

if exist('cprodidx','var')
    if isnumeric(cprodidx)
    cprodidx = num2str(cprodidx,'%02i');
    end
else
    error('Specify product in variable cprodidx');
end

substrates = {};
path = './_StrainBooster/_My_Simulations/';
filesinPath = dir(path);
prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*(?=.xls)'], 'match'); % PXX_ABC
substrates = [substrates{~cellfun(@isempty,substrates)}];
prod = prod(~cellfun(@isempty,prod));
if isempty(substrates)
    substrates = {};
end

for metaboliteC = prod
    metabolite = char(metaboliteC{:});
    cnap = cnapBU;
    for xlsf = [cellstr(metabolite), substrates]
        xlsfile = [char(xlsf) '.xls'];
        if ~exist(xlsfile,'file')
            xlsfile = [char(xlsf) '.xlsx'];
            error('file exists but only as xlsx');
        end
        if ~exist(xlsfile,'file')
            return;
        end
    end
    
    substr = {};
    for i = 1:length(substrates)
        substr = [substr; {[path char(substrates(i)) '.xls']}];
    end

%     maxMCSnumPre = 10;
%     maxMCSnum   = 0;
%     maxMCSsize  = 12;
    if ~exist('maxMCSnumPre','var')
        maxMCSnumPre  = 10;
    end
    if ~exist('maxMCSnum','var')
        maxMCSnum  = 0;
    end
    if ~exist('maxMCSsize','var')
        maxMCSsize  = 12;
    end
    if ~exist('timelimit1','var')
        timelimit1  = 86400;
    end
    if ~exist('timelimit2','var')
        timelimit2  = 86400;
    end
    
    filename=['_StrainBooster/_My_Simulations/Solutions/' cnap.path '-2StepCalc-' metabolite '-' datestr(date,'yyyy-mm-dd')];
    cnap.epsilon = 1e-7;
    preprocess = [];
    sub = [];

    %% Calculate cut sets
    cprintf('blue',[metabolite '\n']);
        % Add new reactions to model from xls
    try
        cnap = CNAaddSpecsAndReacsFromFile(cnap,xlsfile);
    catch
        cprintf([0.8 0.6 0.3],[xlsfile ': no reactions were added to model' char(10)]);
    end
    % Load Cut-Set-Calculation parameters from xls
    [T, t, D, d,notknockable,reacMin,reacMax] = CNAgetMCScalcParamXls( cnap, xlsfile);
    cnap.reacMin = reacMin;
    cnap.reacMax = reacMax;
    inverseCutCount = [];

    [mcs,reacNames,reg, ~, cs] = CNAfindRegMCS2(cnap,T,t,D,d,notknockable,maxMCSnumPre,maxMCSnum,maxMCSsize,[filename '.txt'],preprocess,sub,timelimit1,timelimit2,[],inverseCutCount);

    %save([filename '.mat'],'cnap', 'mcs', 'mcsReadable', 'indexUptakeReacs', 'T', 't', 'D', 'd','notknockable');
    save([filename '.mat'],'cnap', 'mcs', 'inverseCutCount', 'T', 't', 'D', 'd','notknockable');
    cprodidx = ['b' cprodidx];
    if ~exist('noeval','var')
        cprodidx = ['b' cprodidx];
        evalMCSpub;
    end
end