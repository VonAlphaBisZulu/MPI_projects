function [cnap, mcs, T, t, D, d,notknockable] = CNAenumeratecMCSxls( cnap, xls_input_filename , maxMCSnum, maxMCSsize,output_filename, preprocess,sub)
%% Add new reactions to model
try
    cnap = CNAaddSpecsAndReacsFromFile(cnap,xls_input_filename);
catch
    warning([xls_input_filename ': nothing was added to model']);
end
                   
%% read out cMCS-Calculation parameters
TableReads = loadSpecReacXLStoStrArray(xls_input_filename);

% ------------------ identify constraint sheet ---------------------------
% Find "sheet" with parameters for cMCS calculation
[nkorow,nkocol,cMCSSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'notknockable')));

[tRrow,tRcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'targetR'));
[dRrow,dRcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'desiredR'));
[gCrow,gCcol,~]      = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'generalConst'));
[RMinrow,RMinCcol,~] = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmin'));
[RMaxrow,RMaxcol,~]  = find(strcmp(strtrim(TableReads(:,:,cMCSSheet)),'Rmax'));

% ================ Target Reactions ===================
lastrow = tRrow+find(strcmp(strtrim(TableReads(tRrow:end,tRcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = tRrow+find(~strcmp(strtrim(TableReads(tRrow:end,tRcol,cMCSSheet)),''),1,'last')-1;
end
targetReacs = [cellstr(TableReads(tRrow+1:lastrow,tRcol:tRcol+1,cMCSSheet)) num2cell(str2double(TableReads(tRrow+1:lastrow,tRcol+2,cMCSSheet)))];

% ================ Desired Reactions ===================
lastrow = dRrow+find(strcmp(strtrim(TableReads(dRrow:end,dRcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = dRrow+find(~strcmp(strtrim(TableReads(dRrow:end,dRcol,cMCSSheet)),''),1,'last')-1;
end
desiredReacs = [cellstr(TableReads(dRrow+1:lastrow,dRcol:dRcol+1,cMCSSheet)) num2cell(str2double(TableReads(dRrow+1:lastrow,dRcol+2,cMCSSheet)))];

% ================ general Constrains ===================
lastrow = gCrow+find(strcmp(strtrim(TableReads(gCrow:end,gCcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = gCrow+find(~strcmp(strtrim(TableReads(gCrow:end,gCcol,cMCSSheet)),''),1,'last')-1;
end
genericConstrains = [cellstr(TableReads(gCrow+1:lastrow,gCcol:gCcol+1,cMCSSheet)) num2cell(str2double(TableReads(gCrow+1:lastrow,gCcol+2,cMCSSheet)))];
targetReacs  = [targetReacs;  genericConstrains];
desiredReacs = [desiredReacs; genericConstrains];

% ================ Rmin ===================
lastrow = RMinrow+find(strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMinrow+find(~strcmp(strtrim(TableReads(RMinrow:end,RMinCcol,cMCSSheet)),''),1,'last')-1;
end
Rmin = [cellstr(TableReads(RMinrow+1:lastrow,RMinCcol,cMCSSheet)) num2cell(str2double(TableReads(RMinrow+1:lastrow,RMinCcol+1,cMCSSheet)))];

% ================ Rmax ===================
lastrow = RMaxrow+find(strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = RMaxrow+find(~strcmp(strtrim(TableReads(RMaxrow:end,RMaxcol,cMCSSheet)),''),1,'last')-1;
end
Rmax = [cellstr(TableReads(RMaxrow+1:lastrow,RMaxcol,cMCSSheet)) num2cell(str2double(TableReads(RMaxrow+1:lastrow,RMaxcol+1,cMCSSheet)))];

% ================ notknockable ===================
lastrow = nkorow+find(strcmp(strtrim(TableReads(nkorow:end,nkocol,cMCSSheet)),''),1,'first')-2;
if isempty(lastrow)
    lastrow = nkorow+find(~strcmp(strtrim(TableReads(nkorow:end,nkocol,cMCSSheet)),''),1,'last')-1;
end
nokoable = cellstr(TableReads(nkorow+1:lastrow,nkocol,cMCSSheet));
notknockable = [];

%% identify not knockable reactions

for ink = 1:length(nokoable)
    if ~isempty(strcmp(strtrim(cellstr(cnap.reacID)), nokoable(ink)))
        notknockable = [notknockable find(strcmp(strtrim(cellstr(cnap.reacID)), nokoable(ink)))];
    else
        warning(['Reaction ' nokoable(ink) ' could not be found and wasn''t set to ''not knockable''.']);
    end
end

%% set Rmin and Rmax

for irmn = 1:length(Rmin)
    if ~isempty(strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))
        cnap.reacMin((strcmp(strtrim(cellstr(cnap.reacID)), Rmin(irmn,1)))) = cell2mat(Rmin(irmn,2));
    else
        warning(['Reaction ' Rmin(irmn,1) ' could not me found. Boundaries were not set.']);
    end
end

for irmx = 1:length(Rmax)
    if ~isempty(strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))
        cnap.reacMax((strcmp(strtrim(cellstr(cnap.reacID)), Rmax(irmx,1)))) = cell2mat(Rmax(irmx,2));
    else
        warning(['Reaction ' Rmax(irmx,1) ' could not me found. Boundaries were not set.']);
    end
end

%% generate t, T, d, D
% generate t, T
t = zeros(size(targetReacs,1),  1);
T = zeros(size(targetReacs,1),  cnap.numr);

for i = 1:size(targetReacs,1)
    % get sign correctly
    signT = 2* (targetReacs{i,2}=='<') * any(ismember(targetReacs{i,1},'-')) -1; 
                                % + : Tv<=t   
                                % - : TV>=t
                                % also respect sign in definition
    
    % fill T and t
    if any(ismember(targetReacs{i,1},'/'))      % yield
        m = strtrim(strsplit(targetReacs{i,1},'/'));
        m = strrep(m,'-','');
        t(i) = 0;
        T(i, strcmp(cnap.reacID,m(1)) ) =  2* (targetReacs{i,2}=='<') -1;
        T(i, strcmp(cnap.reacID,m(2)) ) =  (2* ~xor(targetReacs{i,2}=='<' , any(ismember(targetReacs{i,1},'-'))) -1) * targetReacs{i,3};
    else                                    % simple constraint
        t(i) = signT * targetReacs{i,3};
        T(i, strcmp(strtrim(cellstr(cnap.reacID)),strrep(targetReacs{i,1},'-','')) ) = signT;
    end
end

% Adapt upper bound of target reactions to 1000
for i = 1:size(T,1)
    tR = find(T(i,:));
    if  any(cnap.reacMax(tR)==0)
        cprintf('green',['one of the target reactions has an upper bound of 0, adapted to 1000' '\n']);
        cnap.reacMax(tR(cnap.reacMax(tR)==0))=1000;
    end
end

% generate d, D
d = zeros(size(desiredReacs,1),  1);
D = zeros(size(desiredReacs,1),  cnap.numr);

for i = 1:size(desiredReacs,1)
    % get sign correctly
    signD = 2* (desiredReacs{i,2}=='<') * any(ismember(desiredReacs{i,1},'-')) -1;
                                % + : Tv<=t   
                                % - : TV>=t
                                % also respect sign in definition
   
    % fill T and t
    if any(ismember(desiredReacs{i,1},'/'))      % yield
        m = strtrim(strsplit(desiredReacs{i,1},'/'));
        m = strrep(m,'-','');
        d(i) = 0;
        D(i, strcmp(cnap.reacID,m(1)) ) =  2* (desiredReacs{i,2}=='<') -1;
        D(i, strcmp(cnap.reacID,m(2)) ) =  (2* ~xor(desiredReacs{i,2}=='<' , any(ismember(desiredReacs{i,1},'-'))) -1) * desiredReacs{i,3};
    else                                    % simple constraint
        d(i) = signD * desiredReacs{i,3};
        D(i, strcmp(strtrim(cellstr(cnap.reacID)),strrep(desiredReacs{i,1},'-','')) ) = signD;
    end
end


% ===== set high boundaries to infinite ==========

cnap.reacMin(cnap.reacMin<-999) = -Inf;
cnap.reacMax(cnap.reacMax>999) = Inf;

%% call MCS Enumerator

if ischar(output_filename) &&  ~isempty(output_filename)
    if ~isempty(preprocess)
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,output_filename,preprocess,sub);
    else
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,output_filename);
    end
else
    if ~isempty(preprocess)
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,preprocess,sub);
    else
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize);
    end
end