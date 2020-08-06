function cnap = CNAremoveSpecsAndReacsDefinedInFile( cnap, spec_reac_set_path )
    %%  Delete species, reactions and a fluxmap that are specified in a xls-file
    % cnap                :  CNA project object
    % spec_reac_set_path  :  [char] relative path from CNA root to xls file
    %
    % ---------------------------------------------------------------------
    % delete species, reactions and a fluxmap that are specified in a xls-file
    %
    % The xls-file contains at least the IDs of the species, reactions and can
    % furthermore supply the name of a map that should be removed from the
    % project. Species, Reactions and the Map-path must be defined on
    % seperate sheets. Those sheets are recognized by the header-names
    % (first line of the sheet).
    %
    % ======== Sheet ======== header ======
    %       Species         spec_id         
    % -------------------------------------
    %    	Reactions       reac_id         % header reac_sbdr if you want
    %                                         to reset boundaries and
    %                                         default rate (put 1 in the
    %                                         column and reaction row)
    % -------------------------------------
    %       Map             map_path       
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 19, 2017

    % read out XLS file
    TableReads = loadSpecReacXLStoStrArray( spec_reac_set_path );
    
    [ridrow,ridcol,reacSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'reac_id')));
    [sidrow,sidcol,specSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'spec_id')));
    [~,mppCol,mapSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'map_path')));

    %% ==========================  Reactions ==============================

    % Find "sheet" with reaction list and predefine struct-array to collect
    % all reactions to remove
    reactions = repmat({'', 2},0,1);
    [~,rsbcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_sbdr'));

    if sum(reacSheet == specSheet)~=0
        error('Invalid input file. Ambiguous definitions found. Make sure that spec_id and reac_id occur in different sheets. No changes have been made to the project');
    elseif (~isempty(sidrow) && sum(sidrow~=1)) || (~isempty(ridrow) && sum(ridrow ~= 1)) || (isempty(ridrow) && isempty(sidrow))
        error('Make sure that spec_id and reac_id are located in the header line');
    else
    end

    % read reaction ids to delete
    lastrow = find(strcmp(strtrim(TableReads(:,ridcol,reacSheet)),''),1,'first')-1;
    % Else take last row with content
    if isempty(lastrow)
        lastrow = find(~strcmp(strtrim(TableReads(:,ridcol,reacSheet)),''),1,'last');
    end
    for row = 2:lastrow
        reacName = TableReads(row,ridcol,reacSheet);
        resetBoundAndDRate = 0;
        if ~isempty(char(TableReads(row,rsbcol,reacSheet)))               % is reaction added or are only boundaries and default rate set
            resetBoundAndDRate = str2double(char(TableReads(row,rsbcol,reacSheet)));
        end
        reactions = [reactions; {reacName} {resetBoundAndDRate}];
    end

    %% ===========================  Species ===============================

    % Find "sheet" with species list and predefine struct-array to collect
    % all species to remove

    species = repmat({''},1,0);

    % read species ids to delete
    lastrow = find(strcmp(strtrim(TableReads(:,sidcol,specSheet)),''),1,'first')-1;
    % Else take last row with content
    if isempty(lastrow)
        lastrow = find(~strcmp(strtrim(TableReads(:,sidcol,specSheet)),''),1,'last');
    end
    for row = 2:lastrow
        species = [species TableReads(row,sidcol,specSheet)];
    end

    %% ============================  Maps =================================

    if cnap.has_gui
        [~,mpncol,~]   = find(strcmp(strtrim(TableReads(:,:,mapSheet)),'map_name'));

        map.map_path = char(TableReads(2,mppCol,mapSheet));

        if ~isempty(TableReads(2,mpncol,mapSheet))
            map.map_name = char(TableReads(2,mpncol,mapSheet));
        else
            [~,map.map_name,~] = fileparts(map.map_path);
        end

        if any(strcmp(cnap.maps(:,1),map.map_name))
            removeMap = 1;
        else
            removeMap = 0;
        end
    end

    %% ============================ Remove ================================
    try
        for reac = [reactions{:,1}]
            % sometimes there is a supplementary space bar added to short ids
            reac_index = find(strcmp(strtrim(cellstr(cnap.reacID)),strtrim(reac)));
            if ~isempty(reac_index)
                % If reaction should not be removed, but only original
                % boundaries and default rate should be restored
                if cell2mat(reactions(strcmp([reactions{:,1}],reac),2))
                    [cnap.reacMin(reac_index), cnap.reacMax(reac_index), cnap.reacDefault(reac_index)] = ...
                    CNAgetGenericReactionData(cnap,reac_index,'reacMin_old','reacMax_old','reacDefault_old');                    
                    cnap = CNAremoveGenericReactionData(cnap,reac_index,'reacMin_old','reacMax_old','reacDefault_old');
                    % If also reaction Box was translocated, relocate it to
                    % original map and position
                    if cell2mat(reactions(strcmp([reactions{:,1}],reac),2))==2
                        delete(cnap.reacBoxes(reac_index,4));
                        cnap.reacBoxes(reac_index,4) = 0;
                        % read out original reac box location from notes
                        % and save to reacBoxes-Field
                        cnap.reacBoxes(reac_index,:) = CNAgetGenericReactionData(cnap,reac_index,'reacBox_old');
                        % remove Layout information from reacNotes
                        cnap = CNAremoveGenericReactionData(cnap,reac_index,'reacBox_old');
%                         cnap.reacNotes(reac_index) = strrep(cnap.reacNotes(reac_index),oldBoxLayout,'');
                        % create new reac box
                        flxMapNo = cnap.reacBoxes(reac_index,5);
                        xpos = cnap.reacBoxes(reac_index,2);
                        ypos = cnap.reacBoxes(reac_index,3);
                        zw=cnap.reacBoxes(reac_index,6);
                        zw1=uicontrol('Style', 'edit','Parent',cnap.figs(flxMapNo,1),'String', '###','Units','normalized','HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',char(reac));
                        set(zw1, 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
                            {'check_for_right_click', 'reaceditmask'}, {'reacenr', reac_index}});
                        cnap.reacBoxes(reac_index,4)=zw1;
                        zoom_single_box(cnap.figs(flxMapNo,:),zw1,xpos,ypos,cnap.reacFontSize(flxMapNo),cnap.reacBoxWidth(flxMapNo),cnap.reacBoxHeight(flxMapNo));

                        if isnan(cnap.reacDefault(reac_index))
                            set(zw1,'String','#');
                        else
                            set(zw1,'String',num2str(cnap.reacDefault(reac_index)));
                        end

                        if(zw==2)       % non-editable
                            set(zw1,'Style', 'text');
                        elseif(zw==3)   % non-visible
                            set(zw1,'Visible','off');
                        end
                    end
                else
                    try
                        cnap = CNAdeleteReaction(cnap, reac_index);
                    catch ME
                        error(['Reaction "',reac,'" could not be deleted',char(10),ME.identifier,': ',ME.message]);
                    end
                end
            end
        end
        for spec = species
            % sometimes there is a supplementary space bar added to short ids
            spec_index = find(strcmp(strtrim(cnap.specID),spec));
            if ~isempty(spec_index)
                try
                    cnap = CNAdeleteSpecies(cnap, spec_index);
                catch
                    warning(['Species "',spec.spec_id,'" could not be removed. Maybe it''s still needed for other reactions']);
                end
            end
        end
        if cnap.has_gui
            if removeMap
                if ~any(cnap.reacBoxes(:,5) == find(strcmp(cnap.maps(:,1),map.map_name)))
                    cnap = CNAremoveMap(cnap,map.map_name);
                else
                    warning(['Map "',map.map_name,'" could not be deleted. There are still reactions on this map.']);
                end
            end
        end
    catch ME
        deleteInvalidTextBoxes(cnap);
        error(['Map, Reactions and Species could not be removed. No changes were made to the project',ME.identifier,': ',ME.message]);
    end
    if cnap.has_gui
        deleteInvalidTextBoxes(cnap);
        cnap = CNAgenerateMap(cnap);
    end
    
end