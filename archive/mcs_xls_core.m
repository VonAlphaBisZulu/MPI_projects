if exist('cprodidx','var')
    if isnumeric(cprodidx)
    cprodidx = num2str(cprodidx,'%02i');
    end
else
    error('Specify product in variable cprodidx');
end

startcna(1)
cnap = CNAloadNetwork({'iJO1366';1},1,1);
cnapBU = cnap;

filepath = './_StrainBooster/_My_Simulations/';
filesinPath = dir(filepath);
prod = regexp({filesinPath.name}, ['^P',cprodidx,'_.*(?=.xls)'], 'match'); % PXX_ABC
prod = prod(~cellfun(@isempty,prod));
subs = regexp({filesinPath.name}, 'S(16)_.*(?=.xls)', 'match'); % S16_ABC
subs = subs(~cellfun(@isempty,subs));

for metaboliteC = prod
    metabolite = char(metaboliteC{:});
    cnap = cnapBU;
    xlsfile = [metabolite '.xls'];
    if ~exist(xlsfile,'file')
        xlsfile = [metabolite '.xlsx'];
    end
    if ~exist(xlsfile,'file')
        return;
    end
    if ~exist('maxMCSsize','var')
        maxMCSsize  = 8;
    end
    
    maxMCSnumPre = 0;
    maxMCSnum   = 1e5;
    timelimit1 = Inf;
    timelimit2 = Inf;
    filename=['_StrainBooster/_My_Simulations/Solutions/' cnap.path '-' metabolite '-' datestr(date,'yyyy-mm-dd')];
    preprocess = [];
    sub = [];

    %% Prepare Model
    % Add new reactions to model from xls
    try
        cnap = CNAaddSpecsAndReacsFromFile(cnap,xlsfile);
    catch
        cprintf([0.8 0.6 0.3],[xlsfile ': no reactions were added to model' char(10)]);
    end
    %% Load Cut-Set-Calculation parameters from xls
    [T, t, D, d,notknockable,reacMin,reacMax] = CNAgetMCScalcParamXls( cnap, xlsfile);
    cnap.reacMin = reacMin;
    cnap.reacMax = reacMax;
    
    %% reduce genome scale Model to core model, also adapt other vectors
    load('/scratch/CNA_SVN/_StrainBooster/Reduced_Model/specAndReac.mat') % infos, which reactions are selected to make up the core network
    otherReacs = cellstr(cnap.reacID(~~sum(T),:));
    [cnap, reac, spec] = buildCoreModel(cnap,spec,reac,otherReacs);
    T = T*sparse(reac,1:length(reac),1);
    D = D*sparse(reac,1:length(reac),1);
    notknockable = find(ismember(reac,notknockable));
    inverseCutCount = [];
    %% Calculate MCS
    [mcs,reacNames,reg, ~, cs] = CNAfindRegMCS2(cnap,T,t,D,d,notknockable,maxMCSnumPre,maxMCSnum,maxMCSsize,[filename '.txt'],preprocess,sub,timelimit1,timelimit2,[],inverseCutCount);
    mcs = CNAregMCSEnumerator2(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,[filename '.txt'],noCompressReacs);

    save([filename '.mat'],'cnap', 'mcs', 'T', 't', 'D', 'd','notknockable');
    evalMCSpub;
end