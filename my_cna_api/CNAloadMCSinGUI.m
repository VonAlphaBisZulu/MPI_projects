function cnap = CNAloadMCSinGUI(cnap,mcs,customValue,lb,ub)
%% sets the entry in all reaction boxes for reactions that are defined in 
%  an mcs vector to 0 or a custom value (such as #)
if nargin >= 3 || cnap.has_gui
    if ~exist('customValue','var')
        customValue = 'gui';
    end
    if ischar(customValue) && strcmp('gui',customValue)
        if nargin > 3
            
        end
        
        cnap= set_local_fields(cnap, {'cutsets', mcs(any(mcs'),:)', 'cutsets_rates_names', cnap.reacID(any(mcs'),:),...
            'num_cutsets', size(mcs,2), 'cutsets_equi_rates', [],'cutsetnr',1});

        cnap.local.cutsets_rates=[];
        for i=1:size(cnap.local.cutsets_rates_names,1)
            zw=mfindstr(cnap.reacID,cnap.local.cutsets_rates_names(i,:));

            if (zw==0)
                warndlg(['Error loading cut sets: Reaction ',deblank(cnap.local.cutsets_rates_names(i,:)),' not found. Loading stopped.'],'Loading cut sets: Error');
                cnap.local.num_cutsets=0;
                cnap.local.cutsets=[];
                return;
            end
            cnap.local.cutsets_rates=[cnap.local.cutsets_rates zw];
        end

        zw=setdiff(1:cnap.numr,cnap.local.cutsets_rates)';
        if(length(zw)>0)
            cnap.local.cutsets_rates_fix=[zw zeros(length(zw),1)];
        else
            cnap.local.cutsets_rates_fix=[];
        end

        cnap= show_cutsets(cnap);
        return
    elseif isnumeric(customValue)
        customValue = num2str(customValue);
    else
        error('invalid value for reaction boxes');
    end
end
if ~any(isnan(mcs) | mcs <0) 
    mcs = -mcs;
end
if cnap.has_gui && any(mcs)
    KOrBoxHandles   = cnap.reacBoxes(mcs ==-1,4);
    noKIrBoxHandles = cnap.reacBoxes(isnan(mcs),4);
    KIrBoxHandles   = cnap.reacBoxes(mcs == 1 ,4);
    if exist('customValue','var')
        set([KOrBoxHandles; KIrBoxHandles; noKIrBoxHandles],'String',customValue);
        set([KOrBoxHandles; KIrBoxHandles; noKIrBoxHandles],'Background',[0.7 0.7 0.7]);
    else
        set(KOrBoxHandles,'String',0);
        set(KOrBoxHandles,'Background',[1 0 0]);
        for i = 1:length(KIrBoxHandles)
            set(KIrBoxHandles,'String','#','Background',[0 1 0]);
        end
        set(noKIrBoxHandles,'String',0,'Background',[1 0.8 0.8])
    end
end
end