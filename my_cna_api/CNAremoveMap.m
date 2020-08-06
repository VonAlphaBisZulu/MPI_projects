function cnap = CNAremoveMap( cnap, map, ignore_box_deletions )
    global cnan;
    cnapBU = cnap;
    %% removes a map (temporarily) from a cna project
    %
    % Usage:    CNAremoveMap( cnap, 'unpleasant_map' )
    %
    % Input:
    %   cnap          :   CNA project object
    %   map           :   [char] map name or (int) map index in cnap.maps of
    %                     the map that should be removed from project
    %
    % Output:
    %   cnap          :   CNA project object
    %
    % Remove a map from a cna project. Reaction boxes that are placed on the
    % map are manually reallocated to other maps/positions through a user
    % interaction dialog
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 16, 2017


    if ~isnumeric(map)
        map = find(strcmp(cnap.maps(:,1),map)); % convert name to index
        if isempty(map)
            error('Map name not found in project');
        end
    elseif map > size(cnap.maps,1)
        error('Map index exceeds numer of elements in projects map list');
    end

    %% =========== reassign all reaction boxes to other maps ==============
    %  Therefore open dialogue box(es)
    if ignore_box_deletions
        reacs_on_map = find(cnap.reacBoxes(:,5) == map);
        cnap.reacBoxes(reacs_on_map,5) = -1;
        delete(cnap.reacBoxes(reacs_on_map,4));
        cnap.reacBoxes(reacs_on_map,4) = 0;
    else
        for reacBoxIndex = find(cnap.reacBoxes(:,5) == map)'
            cancel = 0;

            % open dialog box
            scrsize=get_screen_size();
            cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig=figure('Units','pixels',...
                'NumberTitle','off','Name',['reaction box of "' cnap.reacID(reacBoxIndex,:) '" needs to be relocated'],...
                'Position',[scrsize(3)/2-200 scrsize(4)/2-300 400  200 ],'Menubar','none');
            set(cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig, 'CloseRequestFcn', ['close_local_fig(''', cnap.net_var_name,''',''reacdeffig'');']);
            set(gca,'Visible','off');
            set(gca,'Position',[0,0,1,1]);

            % Initialize drop-down menu for networks
            str = '';
            for j=1:(cnap.nummaps)
                if j ~= map
                    if ~isempty(str)
                        str = [str, '|'];
                    end
                    str=[str,'Map ',num2str(j),': ',cnap.maps{j,1}];
                end
            end

            % Dialog box elements

            text('String','Flux map number',...
                'Position',[30 150],...
                'FontSize',11);
            cnap.local.fmapc = uicontrol('Style','popup',...
                'HorizontalAlignment','left',...
                'Position',[30 140 160 25],...
                'String',str,'Callback',...
                @resetXYloc);

            text('String','X-Position',...
                'Position',[250 150],...
                'FontSize',11);
            cnap.local.xposc=uicontrol('Style','edit',...
                'Position',[250, 140, 80, 25],...
                'FontSize',12,...
                'HorizontalAlignment','left',...
                'Callback',@updatexpos);
            text('String','Y-Position',...
                'Position',[250 80],...
                'FontSize',11);
            cnap.local.yposc=uicontrol('Style','edit',...
                'Position',[250, 70, 80, 25],...
                'FontSize',12,...
                'HorizontalAlignment','left',...
                'Callback',@updateypos);

            cnap.local.getxyf=uicontrol('Style','pushbutton',...
                'String','Get x/y-Pos',...
                'FontSize',12,...
                'Position',[55, 80, 100, 35],'Callback',...
                @getxyloc);

            uicontrol('Style','pushbutton',...
                'String','OK',...
                'FontSize',11,...
                'Position',[60 20 100 35],'Callback',@assignReacToOtherMaploc);
            uicontrol('Style','pushbutton',...
                'String','Cancel',...
                'FontSize',11,...
                'Position',[230 20 100 35],'Callback',@closeDiag);

            activeFigureHandle = {};
            updateFigureHandle();

            % Wait for user confirmation
            uiwait(cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig);
            close(cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig);
            cnap= update_after_change(cnap);
            if cancel
                cnap = cnapBU;
                return;
            end
        end
    end

    %% ================ update reaction-to-map mapping ====================
    % if a map from inbetween the cnap.maps list was removed
    for reacBoxIndex = find(cnap.reacBoxes(:,5) >= map)'
        cnap.reacBoxes(reacBoxIndex,5) = cnap.reacBoxes(reacBoxIndex,5)-1;
    end

    %% ====================== make map disappear ==========================
    updateFigureHandle();
    cnan.open_projects.(cnap.net_var_name).gui.handles.('closeFig') = activeFigureHandle;
    close_local_fig(cnap.net_var_name,'closeFig');

    %% ==================== modify project object =========================

    cnap.maps = [cnap.maps(1:(map-1),:) ; cnap.maps((map+1):end,:)];
    cnap.figs = [cnap.figs(1:(map-1),:)   ; cnap.figs((map+1):end,:)];
    cnap.nummaps = cnap.nummaps-1;
    cnap.reacBoxWidth  = [cnap.reacBoxWidth(1:(map-1)) , cnap.reacBoxWidth((map+1):end)];      % Inherited from first map
    cnap.reacBoxHeight = [cnap.reacBoxHeight(1:(map-1)) , cnap.reacBoxHeight((map+1):end)];
    cnap.specBoxWidth  = [cnap.specBoxWidth(1:(map-1)) , cnap.specBoxWidth((map+1):end)];
    cnap.specBoxHeight = [cnap.specBoxHeight(1:(map-1)) , cnap.specBoxHeight((map+1):end)];
    cnap.reacFontSize  = [cnap.reacFontSize(1:(map-1)) , cnap.reacFontSize((map+1):end)];
    cnap.specFontSize  = [cnap.specFontSize(1:(map-1)) , cnap.specFontSize((map+1):end)];
    
    % make a dialog to suggest saving changes
    saveMapCfgDiag(cnap);

    %% ##################  Supplementary functions ########################
    
        %% ===================== assignReacToOtherMaploc ======================
        function assignReacToOtherMaploc( ~,~ )
            % This function is called before a map is deleted and the
            % reaction boxes that were placed on this map need to be relocated. The
            % function generates a dialog where another map and position can be
            % chosen for the reactions.

            %% ----------- change parameters in cnap.ReacBoxes field --------------
            if ~isempty(cnap.local.xposc.String) && ~isempty(cnap.local.yposc.String)
                % Get mapname from figure name - i.e. "Central Metabolism"
                mapselected = cnap.local.fmapc.String(cnap.local.fmapc.Value,:);
                mapselected = strsplit(mapselected,': ');
                mapselected = strtrim(char(mapselected(end)));

                % Edit parameters of reaction boxes (index of new map, x position, y
                % position)
                cnap.reacBoxes(reacBoxIndex,5) = find(string(cnap.maps) == mapselected);
                cnap.reacBoxes(reacBoxIndex,2) = str2double(cnap.local.xposc.String);
                cnap.reacBoxes(reacBoxIndex,3) = str2double(cnap.local.yposc.String);

                %% ------------------ draw reaction on new map --------------------

                % make sure the cnap-current-figure is also the matlab-current-figure
                if(cnap.figs(cnap.local.fmapc.Value,1)~=gcf)
                    figure(activeFigureHandle);
                end

                % generate and set up textbox graphics element
                zw=cnap.reacBoxes(reacBoxIndex,6);
                zw1=uicontrol('Style', 'edit', 'String', '###','Units','normalized','HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',cnap.reacID(reacBoxIndex,:));
                set(zw1, 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
                    {'check_for_right_click', 'reaceditmask'}, {'reacenr', reacBoxIndex}});
                % register textbox handle in cnap.reacBoxes field
                cnap.reacBoxes(reacBoxIndex,4)=zw1;
                % adjust zoom of textbox
                zoom_single_box(cnap.figs(cnap.local.fmapc.Value,:),zw1,cnap.reacBoxes(reacBoxIndex,2),cnap.reacBoxes(reacBoxIndex,3),cnap.reacFontSize(cnap.local.fmapc.Value),cnap.reacBoxWidth(cnap.local.fmapc.Value),cnap.reacBoxHeight(cnap.local.fmapc.Value));

                % text in box is set to reaction's default value
                if isnan(cnap.reacDefault(reacBoxIndex))
                    set(zw1,'String','#');
                else
                    set(zw1,'String',num2str(cnap.reacDefault(reacBoxIndex)));
                end
                % textbox is made ineditible or invisible if indicated
                if(zw==2)       % non-editable
                    set(reac_index,'Style', 'text');
                elseif(zw==3)   % non-visible
                    set(reac_index,'Visible','off');
                end
            else
                msgbox('please select map and position')
            end
            %% -- delete eventually invalid old text boxes --
            deleteInvalidTextBoxes(cnap);
            uiresume(cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig);
        end
    
        %% ======================== getxyloc ==============================
        function getxyloc(~,~)
            % takes the x and y position that was obtained from the
            % position-choosing routine and puts them into the text fields.
            % Is called, when position on another map has been "clicked"
            
            % gf=gcf;

            updateFigureHandle();

            figure(activeFigureHandle.Number);
            [xpos, ypos]=ginput(1);
            set(cnap.local.xposc,'String',num2str(round(xpos)));
            set(cnap.local.yposc,'String',num2str(round(ypos)));
            figure(activeFigureHandle);
        end
    
        %% ======================= resetXYloc =============================
        function resetXYloc(~,~)
            set(cnap.local.xposc,'String','');
            set(cnap.local.yposc,'String','');
        end

        %% =================== updateFigureHandle =========================
        function updateFigureHandle()
            % set "activeFigureHandle" handle to the map indicated in the
            % drop-down menu.
            % ****** can probably be largely improved ******
            % figure handles usually lie in cnap.figs(:,1) and don't need
            % to be traced back from the figure headline

            % first, get name of map that should be deleted
            mapselected = char(cnap.maps(map,1));
            mapselected = strtrim(char(mapselected));

            % if cnap.local.fmapc is a valid handle, make it the active
            % figure handle
            if isfield(cnap.local,'fmapc')
                if ishandle(cnap.local.fmapc)
                    mapselected = cnap.local.fmapc.String(cnap.local.fmapc.Value,:);
                    mapselected = strsplit(mapselected,': ');
                    mapselected = strtrim(char(mapselected(end)));
                end
            end
            
            % get all matlab figure handles
            handles     = findall(0,'type','figure');
            figureList  = [{get(handles,'Number')}', {get(handles,'Name')}'];

            % find handle of figure whose headline fits wih the selected
            % map
            projfigures = find(~cellfun(@isempty,strfind(figureList{:,2},cnap.net_var_name)))';
            if ~isempty(projfigures)
                handleMap = repmat({'dummy', handles(projfigures(1))},0); % initialize dictionary
                for i = projfigures
                    name = figureList{2}{i};
                    r = strtrim(name((strfind(name, ': ')+1):(strfind(name, ' (')-1)));
                    handleMap = [handleMap ; [r {handles(i)}]];
                end
            else
                error('The project doesnt seem to have any figures');
            end

            % set handle
            activeFigureHandle = handleMap{strcmp(handleMap(:,1),mapselected),2};
        end
    
        %% ========================== others ==============================
        function updatexpos(~,~)
            cnap.local.xposc.Value=str2double(cnap.local.xposc.String);
        end
        function updateypos(~,~)
            cnap.local.yposc.Value=str2double(cnap.local.yposc.String);
        end
        function closeDiag(~,~)
            cancel = 1;
            warning('User interception. No changes have been made to the project');
            uiresume(cnan.open_projects.(cnap.net_var_name).gui.handles.reacdeffig);
        end
end