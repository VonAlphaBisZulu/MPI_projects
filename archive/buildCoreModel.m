% Core model is based on CNA - ECGS - core model
% execute before:
%  1. startcna(1)
%  2. cnap = CNAloadNetwork({'iJO1366';1},1,1);
%  3. load('/scratch/CNA_SVN/_StrainBooster/Reduced_Model/specAndReac.mat')

function [cnapCore, reac, spec] = buildCoreModel(cnap,spec,reac,otherReacs)
% cnap: original model
% spec: species that are kept from the original model
% reac: reactions that are kept from the original model
% otherReacs: at vector of reactions that are not yet in the 'feasible' list, but  should be feasible

%% Copy other model parameters
cnapCore.type = cnap.type;

%% Copy species and reactions

CRandM.specID = cellstr(cnap.specID(spec,:));
CRandM.reacID = cellstr(cnap.reacID(reac,:));

% add supplementary reactions and depending reactions
if nargin >=4
    rix = findStrPos(cnap.reacID,otherReacs);
    [specs, reacs] = findSpecs(cnap,CRandM,rix,0,0);
    reac = sort(unique([reac, reacs, rix]));
    spec = sort(unique([spec, specs]));
end

% iJOcomp.specNotes      = cellstr(strcat('[',iJOcobra.metFormulas(spec),'] <',num2str(iJOcobra.metCharge(spec)),'>'));
% iJOcomp.specNotes      = strrep(iJOcomp.specNotes,' ','');
% iJOcomp.specNotes      = strrep(iJOcomp.specNotes,']<','] <');
cnapCore.specNotes      = cnap.specNotes(spec);

% geneMap                = strrep(iJOcobra.grRules(reac),'and',' & ');
% geneMap                = strrep(geneMap,'or',' | ');
% iJOcomp.reacNotes      = cellstr(strcat(iJOcobra.rxnNames(reac),' {',geneMap,'}'));
cnapCore.reacNotes      = cnap.reacNotes(reac);

cnapCore.stoichMat      = cnap.stoichMat(spec,reac); % Stoichiometric Matrix
cnapCore.specID         = cnap.specID(spec,:);
cnapCore.specLongName   = char(strrep(cellstr(cnap.specLongName(spec,:)),' ','_'));
cnapCore.numis          = length(cnapCore.specNotes);
cnapCore.nums           = length(cnapCore.specNotes);
cnapCore.specExternal   = zeros(1,cnapCore.nums);
cnapCore.specInternal   = 1:cnapCore.nums;
cnapCore.reacID         = cnap.reacID(reac,:);
cnapCore.reacDefault    = cnap.reacDefault(reac);
cnapCore.reacVariance   = cnap.reacVariance(reac);
cnapCore.numr           = length(cnapCore.reacNotes);
cnapCore.reacMin        = cnap.reacMin(reac);
cnapCore.reacMax        = cnap.reacMax(reac);
cnapCore.objFunc        = cnap.objFunc(reac);
cnapCore.macroDefault   = [];
cnapCore.nummac         = 0;
cnapCore.macroComposition = [];
cnapCore.macroSynthboxes= [];
cnapCore.macroID        = [];
cnapCore.macroLongName  = [];
cnapCore.macroBoxes     = [];
cnapCore.mue            = [];
cnapCore.epsilon        = 1e-10;


% reacBoxes
cnapCore.reacBoxes      = cnap.reacBoxes(reac,:);
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('DM_',cnapCore.numr,1),2),5) = 2;
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('DM_',cnapCore.numr,1),2),2) = 196;
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('DM_',cnapCore.numr,1),2),3) = 770-(sum(all(cnapCore.reacID(:,1:3)==repmat('DM_',cnapCore.numr,1),2))-1)*14:14:770;
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('EX_',cnapCore.numr,1),2)&cnapCore.reacBoxes(:,5)~=1,5) = 2;
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('EX_',cnapCore.numr,1),2)&cnapCore.reacBoxes(:,5)~=1,2) = 196;
cnapCore.reacBoxes(all(cnapCore.reacID(:,1:3)==repmat('EX_',cnapCore.numr,1),2)&cnapCore.reacBoxes(:,5)~=1,3) = 1:14:14*sum(all(cnapCore.reacID(:,1:3)==repmat('EX_',cnapCore.numr,1),2)&cnapCore.reacBoxes(:,5)~=1);
otherReacs = find(~(all(cnapCore.reacID(:,1:3)==repmat('EX_',cnapCore.numr,1),2)|all(cnapCore.reacID(:,1:3)==repmat('DM_',cnapCore.numr,1),2)|cnapCore.reacBoxes(:,5)==1));
cnapCore.reacBoxes(otherReacs,5) = ceil((56:length(otherReacs)+55)/220)+1;
cnapCore.reacBoxes(otherReacs,2) = 195+280*floor(mod((56:length(otherReacs)+55)/55,4));
cnapCore.reacBoxes(otherReacs,3) = (mod(56:length(otherReacs)+55,55)+1)*14;

function [specs, reacs] = findSpecs(GS,CP,reacParent,exR,exS)
    if exR == 0
        exR = [];
    end
    if exS == 0
        exS = [];
    end
    [nRspecs,~,~] = find(GS.stoichMat(:,reacParent));
    nRspecs = nRspecs';
    nRspecs = nRspecs(~ismember(nRspecs,exS));
    if ~all(findStrPos(cellstr(CP.specID),cellstr(GS.specID(nRspecs,:)))) && ~isempty(nRspecs)
        % exclude specs that are already in compressed Model
        nRspecs = nRspecs(~findStrPos(cellstr(CP.specID),cellstr(GS.specID(nRspecs,:))));
        % find reactions with species
        [~,nSreacs,~] = find(GS.stoichMat(nRspecs,:));
        if length(nRspecs)>1
            nSreacs = nSreacs';
        end
        % exclude parent reaction
        nSreacs = nSreacs(~ismember(nSreacs,reacParent));
        % exclude other, already considered, reactions
        nSreacs = nSreacs(~ismember(nSreacs,exR));
        % exclude reactions that are already contained
        if ~isempty(nSreacs)
            nSreacs = nSreacs(~contains(cellstr(GS.reacID(nSreacs,:)),cellstr(CP.reacID)));
        end
        % select only reactions that produce spec
        stoich = GS.stoichMat(nRspecs,nSreacs);
        prod = (stoich>0)|((stoich<0)&repmat((GS.reacMin(nSreacs)<0)',length(nRspecs),1));
        nSreacs = nSreacs(any(prod>0,1));
        specs = unique([nRspecs, exS]);
        reacs = unique([nSreacs, reacParent, exR]);
        [specsChild, reacsChild] = findSpecs(GS,CP,nSreacs,reacs,specs);
        specs = [specs, specsChild];
        reacs = [reacs, reacsChild];
    else
        specs = [];
        reacs = [];
    end
end
% finds position of string in list
function indices = findStrPos( str , pattern , opts )
% str       the space that is seached
% pattern   the search keyword or pattern
% opts      options: 0 - normal search ||| 'regex' - regex search
if nargin == 2
    opts = 0;
end
% convert str type to cell
switch class(str)
    case 'string'
        str = cellstr(strtrim(str));
    case 'char'
        str = cellstr(str);
    case 'cell'
        
    otherwise
        error('input 1 of findStrPos doesn''t correct type');
end
% convert pattern type to cell
switch class(pattern)
    case 'string'
        pattern = char(pattern);
    case 'char'
        pattern = cellstr(pattern);
    case 'cell'
        
    otherwise
        error('input 2 of findStrPos isn''t correct type');
end
[rows,cols] = size(pattern);

switch opts
    case 0
        for i = 1:(rows*cols)
            ind                      = find(strcmp(str, pattern(i)));
            indices(1:length(ind),i) = ind;
        end
    case 'regex'
        if opts
            for i = 1:(rows*cols)
                match = regexp(str, pattern, 'match');
                ind = find(~cellfun(@isempty,match));
                indices(1:length(ind),i) = ind;
            end
            
        end
    otherwise
        error('define correct option. 0 for basic search or ''regex''');
end
if (any(size(indices)==0))
    indices = zeros(size(pattern));
%     indices = [];
end
end
end