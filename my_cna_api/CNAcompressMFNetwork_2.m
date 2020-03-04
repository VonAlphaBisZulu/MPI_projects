function [redsmat,irrev,reacidx,metidx,cnapcomp]=CNAcompressMFNetwork_2(cnap,protect_reac,protect_spec,rmCR,rmChoke,rational,blocked_reac,nodisp)
%
% ------------------------------------------------
% CellNetAnalyzer API function 'CNAcompressMFNetwork'
% ------------------------------------------------
% --> loss-free network compression of CNA mass-flow networks
%
% Usage: [redsmat,irrev,reacidx,metidx,cnapcomp]=CNAcompressMFNetwork(cnap,protect_reac,protect_spec,...
%								rmCR,rmChoke,rational,blocked_reac,nodisp)
%
% Based on the stoichiometric matrix and reversibilities of the reactions of a CNA mass-flow project,
% this function compresses the reaction network by
%	(i) removing blocked reactions
%	(ii) lumping enzyme (reaction) subsets
%	(iii) lumping choke points ( = metabolites produced / consumed by one unique reaction and
%             consumed / produced by many reactions)
%	(iv) removing conservation relations (= dependent metabolites).
%
% This type of compression is loss-free meaning, for example, that the compressed system has the
% equivalent set of elementary modes etc.
% Certain reactions or/and metabolites can be protected against compression.
% Compression based on efmtool's rational arithmetic can be used optionally.
%
% Input:
%
%   cnap: (mandatory) is a CellNetAnalyzer (mass-flow) project variable.
%         The function accesses the following fields of cnap (see manual):
%
%      cnap.stoichmat: the stoichiometric matrix of the network
%      cnap.numr: number of reactions (columns in cnap.stoichMat)
%      cnap.numis: number of internal species
%      cnap.mue: index of the biosynthesis reaction; can be empty
%      cnap.macroComposition: matrix defining the stoichiometry of
%        the macromolecules with respect to the metabolites (species); matrix
%        element macroComposition(i,j) stores how much of metabolite i (in mmol)
%        is required to synthesize 1 gram of macromolecule j
%      cnap.specInternal: vector with the indices of the internal species
%      cnap.macroDefault: default concentrations of the macromolecules
%      cnap.reacMin: lower boundaries of reaction rates
%       (if reacMin(i)=0 --> reaction i is considered to be irreversible)
%      cnap.reacMax: upper boundaries of reaction rates
%
% Optional Arguments:
%
%  protect_reac: vector of indices of reactions that should be preserved and not be compressed;
%		 default:[]
%
%  protect_spec: vector of indices of metabolites that should be preserved and not be compressed;
%                default:[]
%
%  rmCR:  	 whether(1) dependent metabolites (=conservation relations) are
%                to be removed or not (0); default: 1
%
%  rmChoke:  	 whether(1) compression of choke points ( = metabolites produced / consumed by
%                one unique reaction and consumed / produced by many reactions) should be conducted
%		 or not (0); default: 1
%
%  rational:  	 whether(1) or not (0) rational arithmetic is to be used for detection of blocked
%                reactions and for compression of enzyme subsets and choke point metabolites
%                (uses efmtool's JAVA compression routine); default: 1
%
%  blocked_reac: indices of reactions that can be deleted before compression; default: [];
%
%  nodisp:	 whether (1) or not (0) display of messages is to be switched off
%
%
% The following results are returned:
%
%  redsmat: 	 the stoichiometric matrix of the reduced system. Only internal metabolites are retained.
%
%  irrev: 	 irreversibilities of the reactions in the reduced system (irrev(i) is 1 (0) if reaction i
%		 of the reduced system is irreversible (reversible)).
%
%  reacidx: 	 reacidx(i,j) gives the stoichiometric coefficient of the original reaction i in the
%                lumped reaction j from the reduced system
%
%  metidx: 	 metidx(i) is the index of the i-th metabolite (row) in redsmat
%		 in the original system (row in smat). Allows mapping of the metabolites
%                in the compressed and in the original system.
%
%  cnapcomp: 	 the compressed network as CNA mass-flow project. Contains also names of the lumped
%                reactions and of the retained metabolites. Note that, in contrast to redsmat,
%		 all external metabolites are retained in the CNA project.
%

% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2019 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de

if nargin<8
    nodisp=0;
end

if nargin<7
    blocked_reac=[];
end
if nargin<6
    rational=1;
end
if nargin<5
    rmChoke=1;
end
if nargin<4
    rmCR=1;
end
if nargin<3
    protect_spec=[];
end

if nargin<2
    protect_reac=[];
end

if(rational==1)
    try
        ch.javasoft.metabolic.compress.CompressionMethod.STANDARD;
    catch
        disp('Cannot access metabolic-efm-all.jar; check the Java classpath.');
        disp('Cannot use rational compression routine; will use standard compression.');
        rational=0;
    end
end

model=CNAgetMFNetwork(cnap);  % just to get the right stoichMat with mue reaction filled

numreac=size(model.stoichMat,2);
usereac=ones(1,numreac);
if(~isempty(blocked_reac))
    zw=intersect(blocked_reac,protect_reac);
    if(~isempty(zw) && ~nodisp)
        disp(['The following reactions were specified as blocked AND protected (those will be considered as blocked and removed):']);
        disp(cnap.reacID(zw,:));
    end
    zw1=zeros(1,numreac);
    zw1(protect_reac)=1;
    zw1(blocked_reac)=[];
    protect_reac=find(zw1);
    usereac(blocked_reac)=0;
end
usereacidx=find(usereac);

extmet=find(model.specExternal);
smatext=model.stoichMat(extmet,:); % here take all reactions; also blocked ones!!
smat=model.stoichMat(model.specInternal,usereacidx);
[meto,reaco]=size(smat);
irrev = (model.reacMin(usereacidx)>=0)';
epsilon=cnap.epsilon;

zw=intersect(protect_spec,extmet);
if(~isempty(zw))
    if(~nodisp)
        disp(['The following metabolites were specified as protected but are external metabolites:']);
        disp(cnap.specID(zw,:));
    end
    protect_spec=setdiff(protect_spec,zw);
end

num_pspec=numel(protect_spec);
num_preac=numel(protect_reac);
aux_reac = [];
metidx=1:size(smat,1);

if(num_pspec)
    if(~nodisp)
        disp(' ');
        disp([num2str(num_pspec),' protected internal metabolites:']);
        disp(cnap.specID(protect_spec,:));
        disp(' ');
    end

    protect_spec_idx=zeros(1,cnap.nums);
    protect_spec_idx(protect_spec)=1;
    protect_spec_idx(extmet)=[];
    protect_spec=find(protect_spec_idx);

    zw=find(sum(smat(:,protect_reac)~=0,2)>=2)'; %% anyway protected species
    if(~isempty(zw))
        if(~nodisp)
            disp('Internal species implicitly protected through protected reactions:');
            disp(cnap.specID(cnap.specInternal(zw),:));
        end
        protect_spec=setdiff(protect_spec,zw);
        num_pspec=numel(protect_spec);
    end
end

if(num_pspec)
    extramat2=zeros(size(smat,1),2*num_pspec);

    for i=1:num_pspec
        extramat2(protect_spec(i),2*i-1)=1;
        extramat2(protect_spec(i),2*i)=-1;
    end
    numr=size(smat,2);
    smat=[smat extramat2];
    idx=numr+1:numr+2*num_pspec;
    aux_reac=[idx];
    irrev=[irrev zeros(1,2*num_pspec)]==1;
end

if(~nodisp && num_preac)
    disp(' ');
    disp([num2str(num_preac),' protected reactions:']);
    disp(cnap.reacID(usereacidx(protect_reac),:));
end

if(rational)
    if(~nodisp)
        disp(' ');
        disp('Using rational compression routines ...');
        disp(' ');
    end
    if(num_preac)
        extramat=smat(:,protect_reac);
        numr=size(smat,2);
        smat=[smat extramat extramat];
        aux_reac=[aux_reac,numr+1:numr+2*num_preac];
        irrev=[irrev zeros(1,2*num_preac)]==1;
    end

    if(rmChoke==0)
        method=1;
    else
        method=0;
    end
    [redsmat,reacidx,mm]= compress_rat_efmtool(smat,~irrev,method,epsilon,rmCR==1);
    if(isempty(redsmat) && ~isempty(reacidx))
        redsmat=zeros(1,size(reacidx,2));
        metidx=metidx(1);
    else
        zw=[];
        for i=1:size(mm,1)
            zw1=find(mm(i,:));
            if(~isempty(zw1))
                zw(end+1)=zw1;
            end
        end
        metidx=metidx(zw);
    end
    irrev=any(reacidx(find(irrev),:),1);

else % not rational
    if rmCR && size(smat,1)>1
        [bm, num_crel]= basic_metabolites(smat,epsilon);
        if(~nodisp && num_crel>0)
            disp(['Found ',num2str(num_crel),' conservation relations and removed dependent metabolites.']);
        end
        smat=smat(bm,:);
        metidx=metidx(bm);
    end

    if(~nodisp)
        disp(' ');
        disp('Using standard (non-rational) compression routines ...');
        disp(' ');
    end
    metidxlast=metidx;

    if islogical(irrev) %A# a logical irrev vector can cause errors during later stages
        zw= zeros(size(irrev));
        zw(irrev)= 1;
        irrev=zw;
    end

    redsmat=smat;
    rr=size(smat,2);
    mm=size(smat,1);
    %metidx=1:mm;
    reacidx=eye(rr,rr);

    cR=[protect_reac, aux_reac];

    ready = 0;
    while ~ready && ~isempty(redsmat)
        ready=1;
        %%% remove non-participating metabolites
        totdeg=sum(abs(redsmat)>epsilon,2);
        todelm=find(totdeg==0);
        redsmat(todelm,:)=[];
        metidx(todelm)=[];
        totdeg(todelm)=[];
        mm=mm-length(todelm);
        if(length(todelm))
            %disp(['0-metabolites: ', num2str(length(todelm))]);
        end

        %%% remove (simple) dead-end metabolites
        todelm=find(totdeg==1);
        if(length(todelm))
            for i=1:length(todelm)
                todelr=find(abs(double(redsmat(todelm(i),:)))>epsilon);
                redsmat(:,todelr)=[];
                reacidx(:,todelr)=[];
                irrev(todelr)=[];
                rr=rr-length(todelr);
            end
            redsmat(todelm,:)=[];
            metidx(todelm)=[];
            totdeg(todelm)=[];
            mm=mm-length(todelm);
            ready=0;
            %disp(['simple dead-ends: ', num2str(length(todelm))]);
        end

        %%% remove ("multiple") dead-end metabolites
        if(ready)
            revdeg=sum(abs(redsmat(:,find(~irrev)))>epsilon,2);
            outdeg=sum(redsmat<-epsilon,2);
            indeg=sum(redsmat>epsilon,2);
            irrmets=find(revdeg==0);
            zw=find(outdeg==0);
            todelm=intersect(irrmets,zw)';
            zw=find(indeg==0);
            todelm=[todelm,intersect(irrmets,zw)'];

            if(length(todelm)>0)
                for i=1:length(todelm)
                    todelr=find(abs(redsmat(todelm(i),:))>epsilon);
                    redsmat(:,todelr)=[];
                    reacidx(:,todelr)=[];
                    irrev(todelr)=[];
                    rr=rr-length(todelr);
                end
                redsmat(todelm,:)=[];
                metidx(todelm)=[];
                mm=mm-length(todelm);
                ready=0;
                %disp(['sink/source with multiple connections: ', num2str(length(todelm))]);
            end
        end

        revmat=-redsmat(:,~irrev);

        %% Lump uniquely produced metabolites (including metabolites with degree 2)
        if(ready && rmChoke)
            %uniprodall=find(indeg+sum(revmat>epsilon,2)==1);
            uniprodall=find((indeg+sum(revmat>epsilon,2)==1) | totdeg==2);  %% consider also all simple chains at this point (totdeg==2)
            %cnap.specID(cnap.specInternal(metidx(uniprodall)),:)

            if ~isempty(uniprodall)
                while(~isempty(uniprodall))
                    uniprod=uniprodall(1);
                    zw= redsmat(uniprod,:);
                    zw3=find(abs(zw)>epsilon);
                    zw2=find(zw>epsilon);
                    zw4=find(zw<-epsilon);
                    %if isempty(cR) || totdeg(uniprod)==2 || (~isempty(zw2) && ~any(ismember(find(reacidx(:,zw2)),cR))) ||...
                    %    (isempty(zw2) && ~any(ismember(find(reacidx(:,find(zw<-epsilon))),cR)))

                    if isempty(cR) || ~any(ismember(find(any(reacidx(:,zw3),2)),cR))
                        break;
                        %new
                        %elseif sum(ismember(find(any(reacidx(:,zw3),2)),cR))==1 && (numel(zw3)==2 || ~any(ismember(find(any(reacidx(:,zw2),2)),cR)))
                        %  break
                        %newend
                    else
                        uniprodall(1)=[];
                    end
                end

                if(length(uniprodall))
                    if(isempty(zw2) || isempty(zw4))   %two reactions pointing in/out of the node; at least one must be reversible
                        zz=find(~irrev(zw3)); %zz must exist
                        zz=zw3(zz(1));
                        zw(zz)=-zw(zz);
                        redsmat(:,zz)=-redsmat(:,zz);
                        reacidx(:,zz)=-reacidx(:,zz);
                        zw2=find(zw>epsilon);
                        zw4=find(zw<-epsilon);
                    end
                    tocomb=zw4;
                    for i=1:length(tocomb)
                        fac=abs(zw(tocomb(i))/zw(zw2));
                        redsmat(:,tocomb(i))=redsmat(:,tocomb(i))+fac*redsmat(:,zw2);
                        reacidx(:,tocomb(i))=reacidx(:,tocomb(i))+fac*reacidx(:,zw2);
                        irrev(tocomb(i))=irrev(tocomb(i)) | irrev(zw2);
                    end
                    redsmat(:,zw2)=[];
                    reacidx(:,zw2)=[];
                    irrev(zw2)=[];
                    redsmat(uniprod,:)=[];
                    metidx(uniprod)=[];
                    ready=0;
                    mm=mm-1;
                    rr=rr-1;
                    %disp(['uniquely produced: ', num2str(length(uniprod))]);
                end
            end
        end

        %% Lump uniquely consumed metabolites
        if(ready && rmChoke)
            %A# is it OK to continue without updating revmat?
            uniconsall=find(outdeg+sum(revmat<-epsilon,2)==1);

            if ~isempty(uniconsall)
                while(~isempty(uniconsall))
                    unicons=uniconsall(1);
                    zw=redsmat(unicons,:);
                    zw2=find(zw<-epsilon);
                    %if totdeg(unicons)==2 || isempty(cR) || (~isempty(zw2) && ~any(ismember(find(reacidx(:,zw2)),cR))) ||...
                    %    (isempty(zw2) && ~any(ismember(find(reacidx(:,find(zw>epsilon))),cR)))
                    zw3=find(abs(zw)>epsilon);
                    if isempty(cR) || ~any(ismember(find(any(reacidx(:,zw3),2)),cR))
                        break;
                        %elseif sum(ismember(find(any(reacidx(:,zw3),2)),cR))==1 && ~any(ismember(find(any(reacidx(:,zw2),2)),cR))
                        %  break;
                    else
                        uniconsall(1)=[];
                    end
                end

                if(length(uniconsall))
                    if(isempty(zw2))
                        zw(find(~irrev))=revmat(unicons,:);
                        zw2=find(zw<-epsilon);
                        redsmat(:,zw2)=-redsmat(:,zw2);
                        reacidx(:,zw2)=-reacidx(:,zw2);
                    end
                    tocomb=find(zw>epsilon);
                    for i=1:length(tocomb)
                        fac=abs(zw(tocomb(i))/zw(zw2));
                        redsmat(:,tocomb(i))=redsmat(:,tocomb(i))+fac*redsmat(:,zw2);
                        reacidx(:,tocomb(i))=reacidx(:,tocomb(i))+fac*reacidx(:,zw2);
                        irrev(tocomb(i))=1;
                    end
                    redsmat(:,zw2)=[];
                    reacidx(:,zw2)=[];
                    irrev(zw2)=[];
                    redsmat(unicons,:)=[];
                    metidx(unicons)=[];
                    mm=mm-1;
                    rr=rr-1;
                    ready=0;
                    %disp(['uniquely consumed: ', num2str(length(unicons))]);
                end
            end
        end

        if ready && ~isempty(redsmat)
            %Compress enzyme subsets
            curprot=find(any(reacidx(cR,:),1));

            [redsmat, sub, irrev, rdind]= subsets_reduction(redsmat,irrev==1, [],curprot,epsilon);

            metidx=metidx(rdind);
            %do null space analysis only once!!
            %if(size(reacidx,2)~=size(sub',2))
            %	ready=0;
            %end
            reacidx=reacidx*(sub');

            %if(~isempty(redsmat) && isempty(metidx))
            %	metidx=metidxlast(1);
            %end
        end

    end %# while(~ready & ~isempty(redsmat))

    if(isempty(redsmat) && ~isempty(reacidx))
        redsmat=zeros(1,rr);
        metidx=metidxlast(1);
    elseif(isempty(redsmat))
        redsmat=[];
        irrev=[];
        metidx=[];
    else
        redsmat(find(abs(redsmat)<epsilon))=0;
    end

    % redsmat= double(redsmat);
    reacidx= double(reacidx);

    if rmCR && size(redsmat,1)>1 %once more; can be useful
        [bm, num_crel]= basic_metabolites(redsmat,epsilon);
        if(~nodisp && num_crel)
            disp(['Found ',num2str(num_crel),' additional conservation relations and removed dependent metabolites.']);
        end
        redsmat=redsmat(bm,:);
        metidx=metidx(bm);
    end
end % if rationale

irrev=double(irrev); %% important for returning irrev

if(~isempty(aux_reac))
    for i=1:numel(aux_reac)
        delr(i)=find(reacidx(aux_reac(i),:));
    end

    redsmat(:,delr)=[];
    irrev(delr)=[];
    reacidx(:,delr)=[];
    reacidx(aux_reac,:)=[];
end

if(~isempty(blocked_reac))  %reinsert
    zw=zeros(numreac,size(reacidx,2));
    zw(usereacidx,:)=reacidx;
    reacidx=zw;
end;


if(nargout<=4)
    if(~nodisp);
        disp(' ');
        disp(['Size of original network: ',num2str(meto),' internal metabolites, ',num2str(reaco),' reactions']);
        disp(['Size of compressed network: ',num2str(size(redsmat,1)),' internal metabolites, ',num2str(size(redsmat,2)),' reactions']);
        disp(' ');
    end
else

    cnapcomp.stoichMat=redsmat;
    cnapcomp.numis=size(cnapcomp.stoichMat,1);
    cnapcomp.nums=cnapcomp.numis+numel(extmet);
    extmetsmatnew=smatext*reacidx;

    cnapcomp.stoichMat=[cnapcomp.stoichMat;extmetsmatnew];
    cnapcomp.specExternal=zeros(1,cnapcomp.nums);
    cnapcomp.specExternal(cnapcomp.numis+1:end)=1;
    cnapcomp.specInternal=find(~cnapcomp.specExternal);
    cnapcomp.specID=[cnap.specID(cnap.specInternal(metidx),:);cnap.specID(extmet,:)];

    reacID=[];
    reacMin=[];
    reacMax=[];

    for i=1:size(reacidx,2)
        zw=find(reacidx(:,i)>0);
        zw2=find(reacidx(:,i)<0);

	lb=[];
	ub=[];

	if(~isempty(zw))
        	lb=max(cnap.reacMin(zw)./reacidx(zw,i));
        	ub=min(cnap.reacMax(zw)./reacidx(zw,i));
	end
	if(~isempty(zw2))
        	lb=max([lb;cnap.reacMax(zw2)./reacidx(zw2,i)]);
        	ub=min([ub;cnap.reacMin(zw2)./reacidx(zw2,i)]);
	end
		
	reacMin(end+1)=lb;
	reacMax(end+1)=ub;

        if(numel(zw)>1)
            str=[];
            for k=1:numel(zw)
                str=[str,deblank(cnap.reacID(zw(k),:)),'*'];
            end
            str=[str,'lumped'];
        else
            str=[deblank(cnap.reacID(zw,:))];
        end
        if(i>1)
            reacID=char(reacID,str);
        else
            reacID=str;
        end
    end

    cnapcomp.type=1;
    cnapcomp.reacID=reacID;
    cnapcomp.reacMin=reacMin'; %P Transposed matrix, else CNAgenerateNetwork throws error
    cnapcomp.reacMax=reacMax'; % 
    cnapcomp.mue=[];

    cnapcomp=CNAgenerateMFNetwork(cnapcomp,nodisp);
    if(~nodisp);
        disp(' ');
        disp(['Size of original network: ',num2str(cnap.numis),' internal metabolites, ',num2str(cnap.numr),' reactions']);
        disp(['Size of compressed network: ',num2str(cnapcomp.numis),' internal metabolites, ',num2str(cnapcomp.numr),' reactions']);
        disp(' ');
    end
end
