clear   targetReacs...
        t ...
        T ...
        tR ...
        desiredReacs ...
        dR ...
        d ...
        D ...
        notknockable ...
        preprocess ...
        maxMCSsize ...
        maxMCSnum ...
        filename


%% choose network and prepare Target and desired matrix
cnap = ECC2comp;    % edit for other networks

%% define constrains, target and desired reactions
% boundaries:

%cnap.reacMin(strcmp(strtrim(cnap.reacID), string('R_ATPM')))=3.15;
%cnap.reacMax(strcmp(strtrim(cnap.reacID), string('R_ATPM')))=3.15;
cnap.reacMax(strcmp(strtrim(cnap.reacID), string('R_GlycUp')))=0;
cnap.reacMax(strcmp(strtrim(cnap.reacID), string('R_SuccUp')))=0;

% genericConstrains
genericConstrains   =  {
                            string('R_GlcUp'),   string('>'), 1e-2;...
                            string('R_ATPM'),   string('>'), 3.15
                       };

% targetReacs:
targetReacs         =  [{    
                            %string('R_FormEx'), string('>'), -1e-4
                            %string('R_AcEx'),   string('>'), 1e-4
                            string('R_LacEx'),  string('>'), 1e-4;
                            %string('R_SuccEx'), string('>'), -1e-4
                            %string('R_EthEx/-R_GlcUp'),  string('>'), 1.4
                            string('R_EthEx'),  string('<'), 5;
                       };genericConstrains];

% desiredReacs:
desiredReacs =         [{    
                            string('R_EthEx'),string('>'),-1e-2;
                       };genericConstrains];

notknockable    = [];
          
maxMCSsize  = 3;
maxMCSnum   = 10000;

filename = '';
preprocess = [];


%% generate t, T, d, D
% generate t, T
t = zeros(size(targetReacs,1),  1);
T = zeros(size(targetReacs,1),  cnap.numr);

for i = 1:size(targetReacs,1)
    % get sign correctly
    signT = 2* (targetReacs{i,2}=='<') * contains(targetReacs{i,1},string('-')) -1; 
                                % + : Tv<=t   
                                % - : TV>=t
                                % also respect sign in definition
    targetReacs{i,1} = strrep(targetReacs{i,1},string('-'),string(''));
    
    % fill T and t
    if contains(targetReacs{i,1},string('/'))      % yield
        m = strtrim(strsplit(targetReacs{i,1},string('/')));
        t(i) = 0;
        T(i, strcmp(cnap.reacID,m(1)) ) =  signT;
        T(i, strcmp(cnap.reacID,m(2)) ) =  signT * targetReacs{i,3};
    else                                    % simple constraint
        t(i) = signT * targetReacs{i,3};
        T(i, strcmp(cnap.reacID,targetReacs{i,1}) ) = signT;
    end
end

clear i m signT

% generate d, D
d = zeros(size(desiredReacs,1),  1);
D = zeros(size(desiredReacs,1),  cnap.numr);

for i = 1:size(desiredReacs,1)
    % get sign correctly
    signD = (2* (desiredReacs{i,2}=='<')-1)* (2*~contains(desiredReacs{i,1},string('-'))-1);
                                % + : Tv<=t   
                                % - : TV>=t
                                % also respect sign in definition
    desiredReacs{i,1} = strrep(desiredReacs{i,1},string('-'),string(''));
    
    % fill T and t
    if contains(desiredReacs{i,1},string('/'))      % yield
        m = strtrim(strsplit(desiredReacs{i,1},string('/')));
        d(i) = 0;
        D(i, strcmp(cnap.reacID,m(1)) ) =  signD;
        D(i, strcmp(cnap.reacID,m(2)) ) =  signD * desiredReacs{i,3};
    else                                    % simple constraint
        d(i) = signD * desiredReacs{i,3};
        D(i, strcmp(cnap.reacID,desiredReacs{i,1}) ) = signD;
    end
end

clear i m signD

%% call MCS Enumerator

filename=['ECC_MCS-',char(datetime('now','Format','yyyy-MM-dd')),'.txt'];

if ischar(filename) &&  ~isempty(filename)
    if ~isempty(preprocess)
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename,preprocess,sub);
    else
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename);
    end
else
    if ~isempty(preprocess)
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,preprocess,sub);
    else
        mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize);
    end
end

clear filename preprocess cnap