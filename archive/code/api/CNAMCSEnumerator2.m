function [mcs] = CNAMCSEnumerator2(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename,preprocess,sub,inverseCutCount)
%
% CellNetAnalyzer API function 'CNAMCSEnumerator'
% -->  Computes Constrained Minimal Cut Sets (cMCSs) in mass-flow networks
%
% Usage:  [mcs] = CNAMCSEnumerator(cnap,T,t,D,d,notknockable,maxMCSnum,maxMCSsize,filename)
% 
% Given a mass-flow project (with or without GUI) and a set of 'undesired' (target) 
% flux vectors (defined by matrix T and vector t) and (optionally) a set of 'desired' 
% flux vectors (defined by matrix D and vector d) Minimal Cut Sets are computed 
% fulfilling the following properties: Knocking out the reactions of an MCS (= setting 
% the corresponding rate to zero) ensures that all (target) flux vectors v obeying
% 
% 	cnap.stoichimat * v = 0
%       reac.Min <= v <= v.reacMax
% 	T*v <= t
% 
% will be blocked (are infeasible) whereas at least one flux vector r fulfilling
% 
% 	cnap.stoichimat * r = 0
%       cnap.reacMin <= r <= cnap.reacMax
% 	D*r <= d
% 
% will be kept functional (if D is empty, no such flux vector must exist).
%
% If D and d are non-empty (i.e., if desired flux vectors have been defined) the MCSs 
% computed represent CONSTRAINED MCSs (cMCSs), otherwise unconstrained MCSs. 
%
% Importantly, cnap.reacMin and cnap.reacMax should be set to -inf / 0 / +inf if the flux
% boundaries are not really known (setting, for example, an (arbitrary) high upper boundary 
% of 10000 instead of +inf can drastically lower the speed of the algorithm as these numerical
% values are explicitely taken into account).
%
% In order to run this function it is necessary that both the MATLAB CPLEX and Java CPLEX 
% interfaces work. If CPLEX is installed under /cluster/apps/cplex-124 the commands for 
% this are (see also startup.m file in CNA's home directory):
%      addpath('/cluster/apps/cplex-124/cplex/matlab/');
%      javaaddpath('/cluster/apps/cplex-124/cplex/lib/cplex.jar');
% Additionally, the MATLAB JVM needs to have the CPLEX shared library on its library path
% which must be set up before (!!) starting MATLAB. For MATLAB versions up to 2010b 
% this can be achieved by adding
%	/cluster/apps/cplex-124/cplex/bin/x86-64_sles10_4.1
% to Matlab's librarypath.txt (or javalibrarypath.txt from version 2012 onwards) configuration file 
%(see also manual).
%
% Inputs: 
%   cnap: (mandatory) is a CellNetAnalyzer (mass-flow) project variable representing
%	  a (metabolic) reaction network. You may easily generate such a structure
%         by the CNAgenerateMFNetwork function.
%
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
%       (if reacMin(i)=0 --> reaction i is irreversible)
%      cnap.reacMax: upper boundaries of reaction rates
%
%   T: the matrix specifying (with vector t) the target flux vectors as given above.  
%      T has Dimension numTargetConst x cnap.numr
%
%   t: the vector specifying (with matrix T) the target flux vectors as given above.  
%      t has Dimension numTargetConst x 1
%
%   D: the matrix specifying (with vector d) the desired flux vectors as given above.  
%      D has Dimension numDesiredConst x cnap.numr. D and d can be empty
%
%   d: the vector specifying (with matrix D) the desired flux vectors as given above.  
%      d has Dimension numDesiredConst x 1. D and d can be empty
%
%   notknockable: row vector with indices of reactions which cannot be cut (knocked-out).
%
%   maxMCSnum: maximal number of MCSs to be computed. This number may sometimes
%      be exceeded slightly. (Just delete the last MCSs if you don't want to have more.)
%      When constrained MCSs are computed with maxMCSnum<inf, the function should be obligatorily
%      called with a filename to ensure proper handling of the upper boundary for the
%      number of cMCSs.
%
%   maxMCSsize: maximum size (cardinality) the MCSs may have. 
%
%   filename (optional): the (c)MCSs can be stored after each iteration (first iteration 
%      calculates (c)MCSs of size 1, second iteration (c)MCSs of size 2 etc.). This will 
%      be a bit slower but ensures that intermediate results are avalable even if the 
%      computation is stopped at a later iteration. If no filename is specified, no results
%      will be stored.
%
% The following result is returned:
%  
%   mcs: the computed (constrained) Minimal Cut Sets. Rows: cMCSs, columns: reactions.
%      A '1' in mcs(i,j) indicates the knockout of reaction j in mcs i.

tic;
mcs=[];
if(nargin<9)
	filename=[];
end
if(nargin<10)
	preprocess=1;
	sub=[];
end

if(preprocess==0 && (nargin<11 || isempty(sub)))
	disp('You must provide the ''sub'' matrix when using the function without preprocessing');
	return;
end
if(nargin<8)
	disp(' ');
	disp('Missing argument, at least 8 arguments required. See help.')
	disp(' ');
	return;
end

if(nargin<9 || isempty(filename))
	dosave=0;
else
	dosave=1;
end

if nargin < 12
    inverseCutCount = [];
end
    
st=initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
stirr= cnap.reacMin >= 0;
[m,q]=size(st);

if(preprocess)
	% FVA to find blocked reactions and to check feasibility of target flux vectors and of desired flux vectors
	disp(' ');
	check_via_fva=1;;

	if(check_via_fva)
		disp('Flux Variability Analysis ...');
		cgp= Cplex();
		cgp.Param.emphasis.numerical.Cur= 1;
		cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
		cgp.Model.A= st;
		cgp.Model.ub= cnap.reacMax;
		cgp.Model.lb= cnap.reacMin;
		cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
		cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
		cgp.Model.obj=zeros(q, 1);
		cgp.Model.sense= 'maximize';
		%res= cgp.solve()
		cgp.DisplayFunc=[];
	
		fva= {cgp.Model.lb, cgp.Model.ub};
		sense= {'minimize', 'maximize'};
		for i= 1:q
		  cgp.Model.obj(:)= 0;
		  cgp.Model.obj(i)= 1;
		  for j= 1:2
		    cgp.Model.sense= sense{j};
		    x= cgp.solve();
		    if x.status ~= 1 && x.status~=2
		      fprintf('status %d %s for %s, %d\n', x.status,x.statusstring, sense{j}, i);
		    end
		    if x.status == 1 
		      fva{j}(i)= x.objval;
		    end
		  end
		end
	
		fvalb= fva{1};
		fvaub= fva{2};
		fva_tol= cgp.Param.simplex.tolerances.optimality.Min;
		fvalb(abs(fvalb) < fva_tol)= 0;
		fvaub(abs(fvaub) < fva_tol)= 0;
		same= fvalb == fvaub;
		blocked_fva= same & (fvalb == 0);
		disp(['Found and removed ',num2str(sum(blocked_fva)),' blocked reactions via FVA!']);
		disp(' ');
	
		%Check feasibility of target flux vectors
		cgp.addRows(-inf(size(T,1),1),T,t);
		x= cgp.solve();
		disp(' ');
		if x.status == 3
		    disp('Target flux vectors infeasible!');
		    fprintf('status %d %s \n', x.status,x.statusstring, sense{j}, i);
		    disp('Exit due to infeasibility of target flux vectors - no cuts required to block target flux vectors!');
		    return;
		else
		    disp('Target flux vectors feasible!');
		end
	
	else
		blocked_fva=[];
	end
	

	extraConst=zeros(0,q);
	extraConstrhs=zeros(0,1);
	zerorow=zeros(1,q);
	for i=1:q
		if(cnap.reacMin(i)~=0 && cnap.reacMin(i)~=-inf)
			extraConst(end+1,:)=zerorow;
			extraConst(end,i)=-1;
			extraConstrhs(end+1,1)=-cnap.reacMin(i);
		end
		if(cnap.reacMax(i)~=inf)
			extraConst(end+1,:)=zerorow;
			extraConst(end,i)=1;
			extraConstrhs(end+1,1)=cnap.reacMax(i);
		end
	end

	if(numel(extraConst>0))
		disp(['Added ',num2str(size(extraConst,1)),' min/max constraints different from 0 or +/- inf!']);
		T=[T;extraConst];
		t=[t;extraConstrhs];
        %P
% 		if(~isempty(D))
% 			D=[D;extraConst];
% 			d=[d;extraConstrhs];
% 		end
	end

	disp([num2str(numel(t)),' constraints for target flux vectors.']);
	disp([num2str(numel(d)),' constraints for desired flux vectors.']);

	%%% reduction 
	Tinvolved=find(any(abs(T)>cnap.epsilon,1));
	Dinvolved=find(any(abs(D)>cnap.epsilon,1));
	TDinvolved=unique([Tinvolved,Dinvolved]);
	notknockable2=unique([notknockable,TDinvolved]);
	
	%%Old version of network compression
	%removal of conservation relations
	%[zw, bc]= rref(st',cnap.epsilon);
	%disp(['Removed ',num2str(size(st,1)-numel(bc)),' metabolites from conservation relations.']);
	%st= st(bc, :);
	%
	%[rd, sub, irrev_rd, sub_irr_viol]= subsets_reduction(st, stirr,blocked_fva, notknockable2);

	%removal of conservation relations
	%[zw, bc]= rref(rd',cnap.epsilon);
	%disp(['Removed ',num2str(size(rd,1)-numel(bc)),' metabolites from conservation relations.']);
	%rd= rd(bc, :);

        %new compression
	%cnap1.stoichMat=st;
  	%cnap1.reacMin=-ones(size(st,2),1);
  	%cnap1.reacMin(find(stirr))=0;
  	%cnap1.reacID=cnap.reacID;
  	%cnap1.epsilon=cnap.epsilon;
  	%cnap1=CNAgenerateMFNetwork(cnap1);
	%[rd, irrev_rd, sub]= CNAcompressMFNetwork(cnap1,notknockable2,[],1,0,1,blocked_fva,0);
	[rd, irrev_rd, sub]= CNAcompressMFNetwork(cnap,notknockable2,[],1,0,1,find(blocked_fva),0);

	sub=sub';  %Transposed version required here (compatible with metatool)

	disp(' ');
	disp(['Final size of reduced system: ',num2str(size(rd))]);
	disp(' ');

	q=size(rd,2);

	todel=[];
	Tidx=[];
	for i=1:numel(Tinvolved)
		zw=find(sub(:,Tinvolved(i)));
		if(isempty(zw))
			%disp(['Warning: reaction ',deblank(cnap.reacID(Tinvolved(i),:)),' contained in definition of target flux vectors is blocked.']);
			%disp(' ');
			todel=[todel, i];
		else
			Tidx=[Tidx zw];
		end
	end
	Tinvolved(todel)=[];
	Tsub=zeros(size(T,1),q);
	Tsub(:,Tidx)=T(:,Tinvolved);

	zw=all(Tsub==0,2);
	todel=[];
	for j=1:numel(zw)
		if(zw(j))
			%if(t(j)==0)
			if(t(j)>=0) %% this happens for instance if a reaction with flux bounds is blocked
				todel=[todel j];
       	        		%disp(['Removed always fulfilled constraint for target flux vectors.']);
				%disp(' ');
			else
       	        		disp(['One constraint for target flux vectors is infeasible. Empty target set! Exit!']);
				disp(' ');
				return;
			end
		end
	end
	Tsub(todel,:)=[];
	t(todel)=[];
	
	if(~isempty(D))
		todel=[];
		Didx=[];
		for i=1:numel(Dinvolved)
			zw=find(sub(:,Dinvolved(i)));
			if(isempty(zw))
       		         	%disp(['Warning: reaction ',deblank(cnap.reacID(Dinvolved(i),:)),' contained in definition of desired flux vectors is blocked.']);
				%disp(' ');
				todel=[todel, i];
       		 	else
				Didx=[Didx zw];
			end
		end
		Dinvolved(todel)=[];
		Dsub=zeros(size(D,1),q);
		Dsub(:,Didx)=D(:,Dinvolved);

		zw=all(Dsub==0,2);
		todel=[];
		for j=1:numel(zw)
			if(zw(j))
				%if(d(j)==0)
				if(d(j)>=0) %% this happens for instance if a reaction with flux bounds is blocked
					todel=[todel j];
       		         		%disp(['Removed always fulfilled constraint for desired flux vectors.']);
					%disp(' ');
				else
       		         		disp(['One constraint for desired flux vectors is infeasible. Exit!']);
					disp(' ');
					return;
				end
			end
		end
		Dsub(todel,:)=[];
		d(todel)=[];
	else
		Dsub=[];
	end

	disp([num2str(numel(t)),' remaining constraints for target flux vectors.']);
	disp([num2str(numel(d)),' remaining constraints for desired flux vectors.']);
	disp(' ');

	notknockablesub=[];
	numk=0;
	for i=1:numel(notknockable)
		zw=find(sub(:,notknockable(i)));
		if(numel(find(sub(zw,:)))==1)
			numk=numk+1;
			notknockablesub(numk)=zw;
		end
	end

	%Check feasibility of desired flux vectors in reduced system
	if(~isempty(Dsub))
       		 disp('Checking feasibility of desired flux vectors ...');
		disp(' ');
		cgp=Cplex();
		cgp.Param.emphasis.numerical.Cur= 1;
		cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
		cgp.Model.A= rd;
		cgp.Model.ub= Inf(size(cgp.Model.A, 2), 1);
		cgp.Model.lb= -Inf(size(cgp.Model.A, 2), 1);
		cgp.Model.lb(irrev_rd ~= 0)= 0;
		cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
		cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
		cgp.Model.obj=zeros(q, 1);
		cgp.Model.sense= 'maximize';
		cgp.DisplayFunc=[];
%		cgp.addRows(-inf(size(Tsub,1),1),-Tsub,-t);  % new!
		cgp.addRows(-inf(size(Dsub,1),1),Dsub,d);
	
		x= cgp.solve();
		disp(' ');
		if x.status == 3 
			    disp('Desired flux vectors infeasible!');
			    fprintf('status %d %s \n', x.status,x.statusstring);
			    disp('Exit due to infeasibility of desired flux vectors!');
    			    return;
		else
			    disp('Desired flux vectors feasible!');
		end
		disp(' ');

		%FVA to find essential reactions for desired flux vectors
		ess_desired=[];
		%fp=fopen('essreacs.txt','w');
		sense= {'minimize', 'maximize'};
	
		for i= 1:size(rd,2)
	  		ess=0;
      		   	cgp.Model.obj(:)= 0;
			cgp.Model.obj(i)= 1;
			val=[NaN NaN];
			for j= 1:2
				cgp.Model.sense= sense{j};
				x= cgp.solve();
      		     		if x.status ~= 1 && x.status~=2 && x.status~=4
     		       			fprintf('Error while testing esential reactions for desired flux vectors.');
		       			fprintf(' Optimization status %d %s for %s, %d\n', x.status,x.statusstring, sense{j}, i);
	       				return;
				end
				if(x.status==4)
					disp(['Reaction ',num2str(i),': Warning: solver status is 4. Assume unbounded solution.']);
		     		elseif(x.status==1 && ((j==1 && x.objval>fva_tol) || (j==2 && x.objval<-fva_tol)))
					ess=j;
					val(j)=x.objval;
				end
     		   	end
		   	if(ess)
				%zw=find(sub(:,i));
				%fprintf(fp,[cnap.reacID(zw,:),' ',num2str(ess),' ',num2str(val(1)),' ',num2str(val(2)),'\n']);
				ess_desired=[ess_desired i];
		   	end
		end

		%fclose(fp);
		disp(['Found ',num2str(numel(unique(ess_desired))),' essential reactions in reduced system for desired flux vectors!']);
%cnap.reacID(ess_desired,:)
 		disp(' ');
		notknockablesub=unique([notknockablesub ess_desired]);
	end

	disp(['Final number of reactions that cannot be knocked-out: ',num2str(numel(notknockablesub))]);
	disp(['Final number of reactions that can be knocked-out: ',num2str(q-numel(notknockablesub))]);
	disp(' ');
	disp('Preprocessing finished!');
	disp(' ');

else %without preprocessing; the T, D etc. must be prepared accordingly; reacMin and reacMax are not considered separately
	Dsub=D;
	d=d;
	Tsub=T;
	t=t;
	rd=initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
	irrev_rd=(cnap.reacMin>=0)';
	q=size(rd,2);
	notknockablesub=notknockable;
	disp(' ');
	disp(['Size of network: ',num2str(size(rd))]);
	disp(['Size of desired system: ',num2str(size(Dsub))]);
	disp(['Size of target system: ',num2str(size(Tsub))]);
	disp(' ');

end


%% MCS computation

if (q-numel(notknockablesub))==0 
	mcs=[];
	disp('No reaction can be knocked-out - no cut set can exist!');
	return;
end


cuts= true(1, q);
cuts(notknockablesub)=false;

%options.cuts= cuts;
%options.method= 5;
%options.workmem= 45056;

%Preparing tests for constrained cut sets
if(~isempty(Dsub))
	cgp=Cplex();
	cgp.Param.emphasis.numerical.Cur= 1;
	cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
  cgp.Param.simplex.tolerances.feasibility.Cur= cgp.Param.simplex.tolerances.feasibility.Min;
	cgp.DisplayFunc=[];
	q= size(rd, 2);
	cgp.Model.A= rd;
	cgp.Model.ub= Inf(q, 1);
	cgp.Model.lb= -Inf(q, 1);
	cgp.Model.lb(irrev_rd ~= 0)= 0;
	cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
	cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);

	cgp.Model.A(end+1:end+size(Dsub,1),:)=Dsub;
	cgp.Model.lhs(end+1:end+size(Dsub,1))=-Inf;
	cgp.Model.rhs(end+1:end+size(Dsub,1))=d;

	cgp.Model.obj=zeros(q, 1);
	cgp.Model.sense= 'maximize';
end
	

cplex_inner= setup_cplex_inner_class_access();
%obj= ConstrainedMinimalCutSetsEnumerator(st, irr, targets, inh, ubi, cuts, flux_lb, flux_ub, des, db)
obj= ConstrainedMinimalCutSetsEnumerator2(rd,irrev_rd,[],Tsub,t,cuts,[],[],[],[],inverseCutCount);
obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 4096);
%obj.cpx.setParam(cplex_inner.DoubleParam.WorkMem, 45056);
if(dosave)
	mcs_cont=logical([]);
	fp=fopen(filename,'w');
	fclose(fp);
	for i=1:maxMCSsize
        obj.ev_size_lb = i-length(inverseCutCount);
		[obj, mcsn, status, comp_time]= shortest_minimal_cut_sets2(obj,i,Inf,Inf,2,[],[],inverseCutCount); %cannot take maxMCSnum here since many will be filtered afterwards
        
        [mcsn,~,~] = unique(mcsn','rows','legacy'); % eliminate redundant mcs
        mcsn=mcsn';

		disp(' ');
		disp([num2str(size(mcsn,2)),' (compressed) MCSs of size ',num2str(i),' found.']);

		%% LP check to filter cMCS from MCS

		if(~isempty(mcsn))
		   if(~isempty(Dsub))
			disp(' ');
			disp('Testing which MCSs are constrained MCSs ...');
			ct= cputime;
			[feasible, objmax]= lp_check_cmcs(mcsn, cgp);
			disp(cputime - ct);
            
			disp(' ');
			disp([num2str(sum(feasible)),' (compressed) constrained MCSs of size ',num2str(i),' found.']);
		
			newmcs=expand_mcs(mcsn(:, feasible), sub)';
			disp(' ');
			disp([num2str(size(newmcs,1)),' (uncompressed) constrained MCSs of size ',num2str(i),' found.']);
			disp(' ');
			if(~isempty(newmcs))
				mcs= [mcs;newmcs];
			end
		   else
			newmcs=expand_mcs(mcsn, sub)';
			disp([num2str(size(newmcs,1)),' (uncompressed) MCSs of size ',num2str(i),' found.']);
			disp(' ');
			if(~isempty(newmcs))
				mcs= [mcs;newmcs];
			end
		   end

		   fp=fopen(filename,'w');
		   fprintf(fp,['MCS up to size ',num2str(i),':\n']);

		   for j=1:size(mcs,1)
			zw=find(mcs(j,:));
			for k=1:numel(zw)
				fprintf(fp,[deblank(cnap.reacID(zw(k),:)),'	']);
			end
			fprintf(fp,'\n');
		   end
		   fclose(fp);
		end
		disp(' ');
		disp([num2str(size(mcs,1)),' MCSs up to size ',num2str(i),' found so far.']);
		disp(' ');
        if(size(mcs,1)>=maxMCSnum)
			break;
        end
        if strcmp(getpref('Internet','SMTP_Server'),'group.mpi-magdeburg.mpg.de') && i>=6
            [pth,fle,~]=fileparts(filename);
            sendmail('schneiderp@mpi-magdeburg.mpg.de',['Otto: ' fle ' cMCS size ' num2str(i) ' completed'],...
            [[num2str(size(mcs,1)),' cMCSs up to size ',num2str(i),' found.'] 10 [num2str(floor(toc/60/60/24)) ' days ' num2str(floor(mod(toc/60/60,24))) 'h' num2str(floor(mod(toc/60,60))) 'm' num2str(floor(mod(toc,60))) 's']]);
            save([pth '/' fle '_cMCS' '.mat'], 'mcs', 'cnap');
        end
	end

%%old	ksmcs.clearModel();
%%old	ksmcs.end();
	disp(' ');
	if(isempty(D))
		disp(['Final result: ',num2str(size(mcs,1)),' MCSs found.']);
	else
		disp(['Final result: ',num2str(size(mcs,1)),' constrained  MCSs found.']);
	end

else
	[obj, mcsn, status, comp_time]= shortest_minimal_cut_sets2(obj,maxMCSsize,maxMCSnum,inf,2,[],[],inverseCutCount);
%old	[mcsn, size_lb]= k_shortest_mcs(maxMCSnum, maxMCSsize, rd, irrev_rd, [], Tsub, t, [], 1, options);

	disp(' ');
	disp([num2str(size(mcs,2)),' (compressed) MCSs found.']);

	%% LP check to filter cMCS from MCS

	if(~isempty(Dsub))
		disp(' ');
		disp('Testing which MCSs are constrained MCSs ...');
		ct= cputime;
        [mcsn,~,~] = unique(mcsn','rows','legacy');
        mcsn = mcsn';
		[feasible, objmax]= lp_check_cmcs(mcsn, cgp);
		disp(cputime - ct);
	
		disp(' ');
		disp([num2str(sum(feasible)),' (compressed) constrained MCSs found.']);

		mcs= expand_mcs(mcsn(:, feasible), sub)';
		disp(' ');
		disp([num2str(size(mcs,1)),' (uncompressed) constrained MCSs found.']);
	else
		mcs= expand_mcs(mcsn, sub)';
		disp([num2str(size(mcs,1)),' (uncompressed) MCSs found.']);
	end
end
toc