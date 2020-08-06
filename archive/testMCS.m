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

%% define target and desired reactions

targetReacs     = { 'R_FormEx',0;...
                    'R_AcEx',0
                  };
          
desiredReacs    = { 'Growth',1e-2;...
                    'R_EthEx',1e-2;...
                    'R_ATPM',3.15;...
                    'R_GlcUp',10
                  };
notknockable    = [];
          
maxMCSsize  = 9;
maxMCSnum   = 5;

filename = '';
preprocess = [];

          
%% generate t, T, d, D
t = cell2mat(targetReacs(:,2));
T = zeros(size(targetReacs,1),  cnap.numr);

for tR = targetReacs(:,1)'
    T(  strcmp(targetReacs(:,1),tR),...
        strcmp(strtrim(cnap.reacID),tR) ) = 1;
end
clear tR targetReacs

d = cell2mat(desiredReacs(:,2));
D = zeros(size(desiredReacs,1),  cnap.numr);

for dR = desiredReacs(:,1)'
    D(  strcmp(desiredReacs(:,1),dR),...
        strcmp(strtrim(cnap.reacID),dR) ) = -1;
end
clear dR desiredReacs


if ischar(filename) &&  ~isempty(filename)
    if ~isempty(preprocess)
        CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename,preprocess,sub)
    else
        CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename)
    end
else
    if ~isempty(preprocess)
        CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,preprocess,sub)
    else
        CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize)
    end
end