function cnap = CNAaddReactionMFNext( cnap, reac, equation, lb, ub, objCoeff, defaultRate, measVar, notes, flxMapNo, editable, xpos, ypos,sbdr)
% CellNetAnalyzer API function 'CNAaddReactionMFN'
% -----------------------------------------------
% Add a reaction to a CNA mass-flow project.
%
% Usage:    cnap = CNAaddReactionMFN( cnap, reac, equation, lb, ub, objCoeff, defaultRate, measVar, notes, flxMapNo, editable, xpos, ypos, sbdr)
% 
% Input: Two input modes:
% 
% 
%   Two input modes:
%
%	--> 1. cnap = CNAaddReactionMFN(cnap, reac)
%   cnap          :   CNA project object
%   [struct] reac: that contains all of the following fields:
%             - [char] reac: reaction identifier (Usually enzyme
%                            abbreviation in capital letters)
%             - [char] equation: i.e. 'A + 2 B = C'
%             - (d)    lb: lower bound for reaction
%             - (d)    ub: upper bound for reaction
%             - (d)    objCoeff
%             - [char] defaultRate
%             - (d)    measVar
%             - [char] notes
%             - (int)  flxMapNo: Index of Flux Map
%             - (int)  editable
%             - (d)      xpos
%             - (int)   sbdr: 1 - set boundaries and default rate only
%                             2 - set boundaries and default rate and place
%                                 relocate reaction box
% 
%   --> 2. CNAaddReactionMFNMFN( cnap, r_id, r_eq, ... (optional params) ) ---
%   reac     :   [char]   Reaction id or name
%   equation :   [char]   Reaction equation (string is also possible)
% 
%      ----- optional parameters ----
%      === name ============= description ====================================== default value =========
%       lb          :   (double) Reaction lower bound                               0
%       ub          :   (double) Reaction upper bound                              100
%       objCoeff    :   (double) Coefficient of the reation in objective function   0
%       defaultRate :   (double) Default reaction rate                              #
%       measVar     :   (double) Measuring variance of reaction rate                0
%       notes       :   [char]   Supplementary notes                                ''
%       flxMapNo    :   (int)    On which Map should the reaction appear         -1 (on generated map)
%                                    (get map counter at: cnap.local.fmapc.Value)
%                                    if -1 is chosen, a new map will be
%                                    generated to place the reaction box upon
%       editable    :   (int)    Textbox is:    1 - editable ;                      1
%                                               2 - non-editable;
%                                               3 - invisible;
%       xpos        :   (double) x position of Textbox on map                       10
%       ypos        :   (double) y position of Textbox on map                       30
%       sbdr        :   (bool)   set boundaries and default value of existing       0 
%                                  reaction instead of adding a new reaction
%                                  0 - new reaction is added
%                                  1 - only parameters are set
%                                  2 - parameters are set and box is relocated
% Output:
% 
%   cnap:         :   CNA project object (with the reaction added)
% 
% 

    %% ======================== handle inputs =============================

    % in case a structure was passed to the function
    if nargin == 2
        equation    =   reac.equation;
        defaultRate =   reac.defaultRate;
        lb          =   reac.lb;
        ub          =   reac.ub;
        objCoeff    =   reac.objCoeff;
        measVar     =   reac.measVar;
        notes       =   reac.notes;
        flxMapNo    =   reac.flxMapNo;
        editable    =   reac.editable;
        xpos        =   reac.xpos;
        ypos        =   reac.ypos;
        sbdr        =   reac.sbdr;
        reac        =   reac.reac;
        % at least the reaction identifier and equation must be defined to add
        % a new reaction.
    elseif nargin < 3
        error([reac ': not enough arguments. At least reaction id and equation must be defined']);
    else
        if nargin < 9
            notes = '';
        end
    end

    % check if default rate is defined
    if (~exist('defaultRate','var'))
        defaultRate='#';
        disp([reac ': No default value for reaction defined! Default rate set to NaN']);
    end

    % convert all variables to valid types if necessary
    reac            =   char(reac);
    equation         =   char(equation);
    defaultRate      =   char(defaultRate);
    notes          =   char(notes);

    % Check and handle inputs, set default values when no parameters are given

    if (~exist('lb','var'))
        lb = 0;
        disp([reac ': No default value for minimal reaction rate defined! Lower bound set to 0']);
    end

    if (~exist('ub','var'))
        ub = 1000;
        disp([reac ': No default value for maximal reaction rate defined! Upper bound set to 1000']);
    end

    if(ub<lb)
        error([reac ': Maximal reaction rate must not be smaller than the minimal rate!']);
    end

    if(~exist('objCoeff','var'))
        disp([reac ': No default value for the coefficient in the objective function defined! Set to 0']);
        objCoeff = 0;
    end

    if(~exist('measVar','var'))
        disp([reac ': No default value for measurement variance defined! Set to 0']);
        measVar = 0;
    end

    if(~exist('notes','var'))
        notes = '';
    end

    if(~exist('flxMapNo','var'))
        disp([reac ': No map defined to place text box! Placed on generated map']);
        flxMapNo = -1;
    end

    if(~exist('editable','var'))
        disp([reac ': Not defined whether textbox is editable! Set to editable']);
        editable = 1;
    end

    if(~exist('xpos','var'))
        disp([reac ': No x-Position of the textbox defined! Set automatically']);
        xpos = 10;
    end

    if(~exist('ypos','var'))
        disp([reac ': No y-Position of the textbox defined! Set automatically']);
        ypos = 30;
    end

    if(~exist('sbdr','var'))
        disp([reac ': Information missing, if reaction is new, or only boundaries and def rate of existing reaction are set']);
    end
    
    %% ================== evaluate reaction equation ======================
    neweq=zeros(1,max(1,cnap.nums));
    cnap.local.noteq=ones(1,max(1,cnap.nums));
    if(strcmp(reac,'mue'))
        if(~isempty(equation))
            error([reac ':Reaction equation for biomass synthesis (mue) is ignored! For defining the stoichiometry of biomass synthesis use instead macromolecule synthesis equations!']);
        end
    else
        if(isempty(equation) && sbdr==0)
            error([reac ':No reaction equation defined!']);
        end
        try
            rem=[equation,' |'];
            [neweq,cnap.local.noteq,~,cnap.local.errval] = read_reaceq(rem,reac,cnap.specID);
            if(~cnap.local.errval)
                return;
            end
        catch ME
            error([reac ': There was an error with the defined reaction equation', char(10),ME.identifier,': ',ME.message]);
        end
    end

    %% ================ add reaction to model object ======================
    % only boundaries and def rate of existing reaction are set
    if exist('sbdr','var') && sbdr
        rindex = find(strcmp(strtrim(cellstr(cnap.reacID)),strtrim(reac)));
        if sbdr == 2 && cnap.has_gui
            delete(cnap.reacBoxes(rindex,4)); % deleting old reaction box
            cnap.reacBoxes(rindex,4) = 0;
            reacBox = sprintf('%.0f,' , cnap.reacBoxes(rindex,:));
            reacBox = reacBox(1:end-1);% strip final comma
            cnap.reacNotes(rindex) = {[char(cnap.reacNotes(rindex)) ' <' reacBox '> ']};
        end
        cnap.reacNotes(rindex) = {[char(cnap.reacNotes(rindex)) '|||' ...
                            'reacMin|' num2str() '||'...
                            'reacMax|' num2str(cnap.reacMax(rindex,1)) '||'...
                            'dRate|'   num2str(cnap.reacDefault(rindex,1))]};
        
        cnap.reacMin(rindex,1)=lb;
        cnap.reacMax(rindex,1)=ub;
        if~(strcmp(defaultRate,'#')==1)
            try
                cnap.reacDefault(rindex,1)=str2double(defaultRate);
                if cnap.has_gui
                    set(cnap.reacBoxes(rindex,4),'String','###');
                    set(cnap.reacBoxes(rindex,4),'Value',str2double(defaultRate));
                end
            catch
                warning(['default rate for "',reac,'" could not be identified. It is set to NaN'])
                cnap.reacDefault(rindex,1)=NaN;
            end
        else
            cnap.reacDefault(rindex,1)=NaN;
            if cnap.has_gui && cnap.reacBoxes(rindex,4) ~= 0
                set(cnap.reacBoxes(rindex,4),'String','#');
                set(cnap.reacBoxes(rindex,4),'Value',0);
            end
        end
        cnap.unsaved_changes= true;
        if sbdr == 1
            return;
        end
    else
        % Abort if not only boundaries should be set, but reaction already
        % exists
        if(mfindstr(cnap.reacID,reac))
            error('The reaction identifier exists already! Use another one!');
        end

        cnap.numr=cnap.numr+1;
        if(cnap.numr==1)
            cnap.reacID=char(reac);
        else
            cnap.reacID=char(cnap.reacID,reac);
        end
        % evaluate default reaction rate
        if~(strcmp(defaultRate,'#')==1)
            try
                cnap.reacDefault(cnap.numr,1)=str2double(defaultRate);
            catch
                warning(['default rate for "',reac,'" could not be identified. It is set to NaN'])
                cnap.reacDefault(cnap.numr,1)=NaN;
            end
        else
            cnap.reacDefault(cnap.numr,1)=NaN;
        end
        % set all reaction parameters
        cnap.objFunc(cnap.numr,1)=objCoeff;
        cnap.reacMin(cnap.numr,1)=lb;
        cnap.reacMax(cnap.numr,1)=ub;
        cnap.reacVariance(cnap.numr,1)=measVar;
        cnap.stoichMat(:,cnap.numr)=neweq';
        cnap.reacNotes{cnap.numr}=str2notes(notes);
        rindex = cnap.numr;
    end

    %% =================== add reaction to gui ============================
    % Section was copied from reacnewmaskeval, line 110ff (Jun 14, 2017)
    % and adapted

    % check if the indicated map-index is in the range of the maps
    % contained in the project. Else generate a map to put the reaction box
    % on.
    
    if isfield(cnap,'reacBoxes')
        % define Reaction box
        if cnap.has_gui
            cnap.reacBoxes(rindex,:)=[rindex xpos ypos 0 flxMapNo editable];
            if flxMapNo <= size(cnap.figs,1) && flxMapNo >= 1

                zw=cnap.reacBoxes(rindex,6);
                zw1=uicontrol('Style', 'edit','Parent',cnap.figs(flxMapNo,1),'String', '###','Units','normalized','HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',reac);
                set(zw1, 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
                    {'check_for_right_click', 'reaceditmask'}, {'reacenr', rindex}});
                cnap.reacBoxes(rindex,4)=zw1;
                zoom_single_box(cnap.figs(flxMapNo,:),zw1,xpos,ypos,cnap.reacFontSize(flxMapNo),cnap.reacBoxWidth(flxMapNo),cnap.reacBoxHeight(flxMapNo));

                if isnan(cnap.reacDefault(rindex))
                    set(zw1,'String','#');
                else
                    set(zw1,'String',num2str(cnap.reacDefault(rindex)));
                end

                if(zw==2)       % non-editable
                    set(zw1,'Style', 'text');
                elseif(zw==3)   % non-visible
                    set(zw1,'Visible','off');
                end
            else
                % generate map if necessary (only when single parameters were
                % passed, so that the number of re-generation stays considerable)
                if nargin > 2
                    cnap = CNAgenerateMap(cnap);
                end
            end
        else
            cnap.reacBoxes(cnap.numr,:)=[cnap.numr xpos ypos -1 -1 editable];
        end
        if(strcmp(reac,'mue'))
            cnap.mue=cnap.numr;
        end
    
    cnap.unsaved_changes= true;

    end
end