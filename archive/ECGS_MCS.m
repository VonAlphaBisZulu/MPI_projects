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
cnap = ECGS;    % edit for other networks

%% define constrains, target and desired reactions
% boundaries:

% finite bounds
cnap.reacMin(cnap.reacMin==-inf)=-10000;
cnap.reacMax(cnap.reacMax==inf)=10000;

%cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_ATPM'))=3.15;
%cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_ATPM'))=3.15;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_Ec_biomass_iJO1366_core_53p95M'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_o2_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_etoh_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_lac_DASH_D_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_ac_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_for_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_succ_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_h2o_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_h2_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_h_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_co2_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_DM_4CRSOL'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_DM_5DRIB'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_DM_AMOB'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_DM_MTHTHF'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_succ_LPAREN_e_RPAREN_'))=10000;
cnap.reacMax(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_o2_LPAREN_e_RPAREN_'))=0;

% glucose limit to -20 and atp limit to 3.15:
cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_o2_LPAREN_e_RPAREN_'))=0;
cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_ac_LPAREN_e_RPAREN_'))=0;
cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_EX_glc_LPAREN_e_RPAREN_'))=-20;
cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_ATPM'))=3.15;
cnap.reacMin(strcmp(strtrim(cellstr(cnap.reacID)), 'R_NADH16pp'))=-10000;
cnap.reacMax(cnap.reacMin(1:332)<0)= 10000;

% genericConstrains
genericConstrains   =  {
                            'R_EX_glc_LPAREN_e_RPAREN_',   '>', -18;...
                            'R_ATPM',   '>', 3.15
                       };

% targetReacs:
targetReacs         =  [{    
                            %'R_FormEx', '>', -1e-4
                            %'R_AcEx',   '>', 1e-4
                            %'R_LacEx',  '>', 1e-4;
                            %'R_SuccEx', '>', -1e-4
                            'R_EX_etoh_LPAREN_e_RPAREN_/-R_EX_glc_LPAREN_e_RPAREN_',  '<', 1.4
                            %'R_EthEx',  '<', 5;
                       };genericConstrains];

% desiredReacs:
desiredReacs =         [{    
                            'R_Ec_biomass_iJO1366_core_53p95M','>',0.05;...
                            'R_EX_etoh_LPAREN_e_RPAREN_/-R_EX_glc_LPAREN_e_RPAREN_',  '>', 1.4
                       };genericConstrains];
          
maxMCSsize  = 3;
maxMCSnum   = 1000000;

preprocess = [];

%% Not knockable Reactions
% Don't knock R_EX-Reactions, don't knock "tex" reactions, don't knock
% diffusion-controlled reactions

% ATP-dependend- and iron-transporters stay knockable
knockOthers = (strcmp(strtrim(cellstr(cnap.reacID)), 'R_ARBTNtex')|...
    strcmp(strtrim(cellstr(cnap.reacID)), 'R_FE3HOXUtex')|...
    strcmp(strtrim(cellstr(cnap.reacID)), 'R_FECRMUtex')|...
    strcmp(strtrim(cellstr(cnap.reacID)), 'R_FEENTERtex')|...
    strcmp(strtrim(cellstr(cnap.reacID)), 'R_FEOXAMUtex'))';
% these lines equal the "contains" function
knockabc=ismember(cnap.reacID','abctex');
knockabc=sum((      knockabc(1:end-5,:)+...
                    knockabc(2:end-4,:)+...
                    knockabc(3:end-3,:)+...
                    knockabc(4:end-2,:)+...
                    knockabc(5:end-1,:)+...
                    knockabc(6:end,:))==6);

% p -> e and e -> ex are not knockable when diffusion-controlled
notKnocktex=ismember(cnap.reacID','tex');
notKnocktex=sum((   notKnocktex(1:end-2,:)+...
                    notKnocktex(2:end-1,:)+...
                    notKnocktex(3:end,:))==3);

notKnockREX=ismember(cnap.reacID','R_EX_');
notKnockREX=sum((   notKnockREX(1:end-4,:)+...
                    notKnockREX(2:end-3,:)+...
                    notKnockREX(3:end-2,:)+...
                    notKnockREX(4:end-1,:)+...
                    notKnockREX(5:end,:))==5);

notKnocktipp=ismember(cnap.reacID','tipp');
notKnocktipp=sum((  notKnocktipp(1:end-3,:)+...
                    notKnocktipp(2:end-2,:)+...
                    notKnocktipp(3:end-1,:)+...
                    notKnocktipp(4:end,:))==4);

% c -> p are not knockable when diffusion-controlled (not always true, but speeds up calculation)

M_h_c=ismember(cnap.specID','M_h_c');
M_h_c=sum((         M_h_c(1:end-4,:)+...
                    M_h_c(2:end-3,:)+...
                    M_h_c(3:end-2,:)+...
                    M_h_c(4:end-1,:)+...
                    M_h_c(5:end,:))==5);
                
M_succ=ismember(cnap.specID','M_succ');
M_succ=sum((        M_succ(1:end-5,:)+...
                    M_succ(2:end-4,:)+...
                    M_succ(3:end-3,:)+...
                    M_succ(4:end-2,:)+...
                    M_succ(5:end-1,:)+...
                    M_succ(6:end,:))==6);
                
hDependend = sum(cnap.stoichMat(M_h_c|M_succ,:)~=0);

notKnocktpp=ismember(cnap.reacID','tpp'); % IMPORTANT! Search is case sensitive
notKnocktpp=sum((   notKnocktpp(1:end-2,:)+...
                    notKnocktpp(2:end-1,:)+...
                    notKnocktpp(3:end,:))==3);

notknockable    = find( notKnocktex&(~knockOthers)&(~knockabc)|...
                        notKnockREX|...
                        notKnocktipp|...
                        notKnocktpp&(~hDependend));
                    
clear M_h_c M_succ notKnocktex knockOthers knockabc notKnockREX notKnocktipp notKnocktpp hDependend

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
    targetReacs{i,1} = strrep(targetReacs{i,1},'-','');
    
    % fill T and t
    if any(ismember(targetReacs{i,1},'/'))      % yield
        m = strtrim(strsplit(targetReacs{i,1},'/'));
        t(i) = 0;
        T(i, strcmp(cnap.reacID,m(1)) ) =  signT;
        T(i, strcmp(cnap.reacID,m(2)) ) =  signT * targetReacs{i,3};
    else                                    % simple constraint
        t(i) = signT * targetReacs{i,3};
        T(i, strcmp(strtrim(cellstr(cnap.reacID)),targetReacs{i,1}) ) = signT;
    end
end

clear i m signT

% generate d, D
d = zeros(size(desiredReacs,1),  1);
D = zeros(size(desiredReacs,1),  cnap.numr);

for i = 1:size(desiredReacs,1)
    % get sign correctly
    signD = (2* (desiredReacs{i,2}=='<')-1)* (2*~any(ismember(desiredReacs{i,1},'-'))-1);
                                % + : Tv<=t   
                                % - : TV>=t
                                % also respect sign in definition
    desiredReacs{i,1} = strrep(desiredReacs{i,1},'-','');
    
    % fill T and t
    if any(ismember(desiredReacs{i,1},'/'))      % yield
        m = strtrim(strsplit(desiredReacs{i,1},'/'));
        d(i) = 0;
        D(i, strcmp(cnap.reacID,m(1)) ) =  signD;
        D(i, strcmp(cnap.reacID,m(2)) ) =  signD * desiredReacs{i,3};
    else                                    % simple constraint
        d(i) = signD * desiredReacs{i,3};
        D(i, strcmp(strtrim(cellstr(cnap.reacID)),desiredReacs{i,1}) ) = signD;
    end
end

clear i m signD

%% call MCS Enumerator

filename=['ECGS_MCS-',datestr(date,'yyyy-mm-dd'),'.txt'];

% if ischar(filename) &&  ~isempty(filename)
%     if ~isempty(preprocess)
%         mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename,preprocess,sub);
%     else
%         mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename);
%     end
% else
%     if ~isempty(preprocess)
%         mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,preprocess,sub);
%     else
%         mcs = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize);
%     end
% end

clear filename preprocess