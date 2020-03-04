function cnap = CNAaddSpecsAndReacsFromFile( cnap, spec_reac_set_path, bIgnoreFalseEntries )
global cnan;
    %%  Import species, reactions and a fluxmap from a xls-file
    % cnap                :  CNA project object
    % spec_reac_set_path  :  [char] relative path from CNA root to xls file
    %
    % ---------------------------------------------------------------------
    % import species, reactions and a fluxmap from a xls-file
    %
    % The xls-file can contain information about species, reactions and can
    % furthermore supply the filepath of a map that should be added to the
    % project. Species, Reactions and the Map-path must be defined on
    % seperate sheets. Those sheets are recognized by the header-names
    % (first line of the sheet). Some header-names are required others are
    % optional. All possible fields are listed below.
    %
    % ======== Sheet ======== header ====== required or optional === description ====
    %       Species         spec_id         REQ
    %       Species         spec_bool_ext   REQ         is Spec. external 0/1
    %       Species         spec_name       *opt
    %       Species         spec_notes      *opt
    % ---------------------------------------------------------------------
    %    	Reactions       reac_id         REQ
    %    	Reactions       reac_eq         REQ
    %    	Reactions       reac_dRate      *opt        default rate
    %    	Reactions       reac_lb         *opt        lower bound
    %    	Reactions       reac_ub         *opt        upper bound
    %    	Reactions       reac_objcoeff   *opt        coeff in objective function
    %    	Reactions       reac_measVar    *opt        measuring variance
    %    	Reactions       reac_notes      *opt
    %    	Reactions       reac_flxMpNo    *opt        index of flux map
    %                                                     to appear on:
    %                                                   1 -> n place on existing map
    %                                                     0    place on map defined in xls-file
    %                                                    -1    place on (newly) generated map
    %    	Reactions       reac_editable   *opt        is reaction box:
    %                                                   editable (1) ineditable (2) or invisible (3)
    %    	Reactions       reac_xpos       *opt        position on map
    %    	Reactions       reac_ypos       *opt
    %       Reactions       reac_sbdr       *opt        if 1, only the default rate and boundaries
    %                                                   of an eventually existing reaction are set
    % ---------------------------------------------------------------------
    %       Map             map_path        REQ         relative path to
    %                                                   fluxmap-imagefile (from CNA root)
    %       Map             map_name        *opt
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 15, 2017

    %% ======= read out xls-file and pre-process information ==============
    % read out xls-file and put data into a 3D-string-array
    TableReads = loadSpecReacXLStoStrArray(spec_reac_set_path);

    % Identify sheets that contain Species, Reaction and Map information

    % ------------------ identify Species sheet ---------------------------
    % Find "sheet" with species list and predefine empty struct-array to collect
    % all species to add
    [sidrow,sidcol,specSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'spec_id')));
    species = repmat(struct(...
                            'spec_id',          '',...
                            'spec_ext',    0,...
                            'spec_name',        '',...
                            'spec_notes',''),0);

    % ------------------ identify Reacions sheet --------------------------
    % Find "sheet" with reaction list and predefine empty struct-array to collect
    % all reactions to add
    [ridrow,ridcol,reacSheet]  = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'reac_id')));
    reactions = repmat(struct(...
                                'reac',           '',...
                                'equation',     '',...
                                'lb',           0,...
                                'ub',           100,...
                                'objCoeff',     0,...
                                'defaultRate',  '#',...
                                'measVar',      0,...
                                'notes',       '',...
                                'flxMapNo',     1,...
                                'editable',     1,...
                                'xpos',         20,...
                                'ypos',         20,...
                                'sbdr',         0),0);

    % ---------------------- identify Map sheet ---------------------------
    % Find "sheet" with path to map and set path
    [~,mppCol,mapSheet]  = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'map_path')));
    map = repmat(struct('map_path','','map_name',''),1);

    if sum(reacSheet == specSheet)~=0
        error('Invalid input file. Ambiguous definitions found. Make sure that spec_id and reac_id occur in different sheets. No changes have been made to the project');
    elseif (~isempty(sidrow) && sum(sidrow~=1)) || (~isempty(ridrow) && sum(ridrow ~= 1)) || (isempty(ridrow) && isempty(sidrow))
        error('Make sure that spec_id and reac_id are located in the header line');
    end

    %% = gather information on map, species and reaction in struct arrays =

    %  ---------------------- Parse map -----------------------------------
    if ~isempty(mapSheet)
        % find column where map name defined
        [~,mpncol,~]   = find(strcmp(strtrim(TableReads(:,:,mapSheet)),'map_name'));
        map.map_path = char(TableReads(2,mppCol,mapSheet));
        [fpath,file,ext] = fileparts(which(spec_reac_set_path));
        map.map_path = [relPathToFrom(fpath,[pwd '/' cnap.path]) '/' map.map_path];

        % set map name (if not defined in xls-file, the map is named after the image file)
        if ~isempty(TableReads(2,mpncol,mapSheet))
            map.map_name = char(TableReads(2,mpncol,mapSheet));
        else
            [~,map.map_name,~] = fileparts(map.map_path);
        end
    end

    %  -------------------- Parse species ---------------------------------
    if ~isempty(specSheet)
        % find columns where species parameters are defined
        %    spec_bool_ext is mandatory
        [sexrow,sexcol,~]   = find(strcmp(strtrim(TableReads(:,:,specSheet)),'spec_bool_ext'));
        
        if isempty(sexrow) || sum(sexrow~=1)
            error('Make sure spec_bool_ext located in the header line');
        end
        [~,snmcol,~]   = find(strcmp(strtrim(TableReads(:,:,specSheet)),'spec_name'));
        [~,sntcol,~]   = find(strcmp(strtrim(TableReads(:,:,specSheet)),'spec_notes'));

        % Read species information, define them and structs and list them in an array
        % Iterate through all ids. Stop at first occurence of empty spec_id field;
        lastrow = find(strcmp(strtrim(TableReads(:,sidcol,specSheet)),''),1,'first')-1;
        % Else take last row with content
        if isempty(lastrow)
            lastrow = find(~strcmp(strtrim(TableReads(:,sidcol,specSheet)),''),1,'last');
        end
        % Iteration through species/lines
        for row = 2:lastrow
            % create dummy structure
            newspec = repmat(struct(...
                'spec_id','',...
                'spec_ext',0,...
                'spec_name','',...
                'spec_notes',''),1);   % Initialize new species
            % fill structure
            newspec.spec_id = TableReads(row,sidcol,specSheet); % ID
            try                                                 % External or internal
                newspec.spec_ext = str2double(char(TableReads(row,sexcol,specSheet)));
            catch ME
                error([newspec.spec_id ': Ambiguous entry in line in xls-file. Aborting. No changes have been made to the project.',char(10),ME.identifier,': ',ME.message]);
            end
            if ~isempty(TableReads(row,snmcol,specSheet))       % species Name
                newspec.spec_name = TableReads(row,snmcol,specSheet);
            else
                warning([char(newspec.spec_id) ': No name defined - metabolie was named by its id']);
                newspec.spec_name = newspec.spec_id;
            end
            if ~isempty(TableReads(row,sntcol,specSheet))       % species Notes
                newspec.spec_notes = TableReads(row,sntcol,specSheet);
            end
            species = [species newspec];
        end
    end
    %  ------------------------ Parse Reactions ---------------------------
    if ~isempty(reacSheet)
        % find columns where reaction parameters are defined
        %    reac_eq and reac_dRate are mandatory
        [reqrow,reqcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_eq'));

        [~,rdrcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_dRate'));
        [~,rlbcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_lb'));
        [~,rubcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_ub'));
        [~,roccol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_objcoeff'));
        [~,rmvcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_measVar'));
        [~,rntcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_notes'));
        [~,rfmcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_flxMpNo'));
        [~,redcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_editable'));
        [~,rxpcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_xpos'));
        [~,rypcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_ypos'));
        [~,rsbcol,~]   = find(strcmp(strtrim(TableReads(:,:,reacSheet)),'reac_sbdr'));

        % Read Reaction information, define them and structs and list them in an array
        % Iterate through all ids. Stop at first occurence of empty reac_id field;
        lastrow = find(strcmp(strtrim(TableReads(:,ridcol,reacSheet)),''),1,'first')-1;
        % Else take last row with content
        if isempty(lastrow)
            lastrow = find(~strcmp(strtrim(TableReads(:,ridcol,reacSheet)),''),1,'last');
        end
        for row = 2:lastrow
            newreac = repmat(struct(...
                'reac',           '',...
                'equation',     '',...
                'lb',           0,...
                'ub',           100,...
                'objCoeff',     0,...
                'defaultRate',  '#',...
                'measVar',      0,...
                'notes',       '',...
                'flxMapNo',     '-1',...
                'editable',     1,...
                'xpos',         20,...
                'ypos',         20,...
                'sbdr',         0),1);  % Initialize new Reaction
            
            newreac.reac = strtrim(char(TableReads(row,ridcol,reacSheet)));         % ID
            try
                newreac.equation  = char(TableReads(row,reqcol,reacSheet)); % Reaction equation
            catch ME
                error([newspec.reac ': Ambiguous entry in line in xls-file. Aborting. No changes have been made to the project.',char(10), ME.identifier,': ',ME.message]);
            end
            if ~isempty(char(TableReads(row,rlbcol,reacSheet)))         % Reaction lower bound
                newreac.lb = str2double(char(TableReads(row,rlbcol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,rubcol,reacSheet)))         % Reaction upper bound
                newreac.ub = str2double(char(TableReads(row,rubcol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,roccol,reacSheet)))         % Objective coefficient
                newreac.objCoeff = str2double(char(TableReads(row,roccol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,rdrcol,reacSheet)))         % Default reaction rate
                newreac.defaultRate = char(TableReads(row,rdrcol,reacSheet));
            end
            if ~isempty(char(TableReads(row,rmvcol,reacSheet)))         % Measuring variance
                newreac.measVar = str2double(char(TableReads(row,rmvcol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,rntcol,reacSheet)))         % Notes
                newreac.notes = char(TableReads(row,rntcol,reacSheet));
            end
            if ~isempty(char(TableReads(row,rfmcol,reacSheet)))         % Flux map index
                newreac.flxMapNo = char(TableReads(row,rfmcol,reacSheet));
            end
            if ~isempty(char(TableReads(row,redcol,reacSheet)))         % Is editable
                if any(str2double(char(TableReads(row,redcol,reacSheet))) == [1,2,3])
                    newreac.editable = str2double(char(TableReads(row,redcol,reacSheet)));
                else
                    newreac.editable = 1;
                end
            end
            if ~isempty(char(TableReads(row,rxpcol,reacSheet)))        % x position
                newreac.xpos = str2double(char(TableReads(row,rxpcol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,rypcol,reacSheet)))        % y position
                newreac.ypos = str2double(char(TableReads(row,rypcol,reacSheet)));
            end
            if ~isempty(char(TableReads(row,rsbcol,reacSheet)))        % is reaction added or are only boundaries and default rate set
                newreac.sbdr = str2double(char(TableReads(row,rsbcol,reacSheet)));
            else
                newreac.sbdr = 0;
            end

            reactions = [reactions newreac];
            if isempty(reqrow)
                if ~newreac.sbdr
                    error('It was attempted to add a new reaction to the model without defining a reaction equation.');
                end
            end
            if prod([double(~isempty(reqrow)) double(reqrow~=1)]) || ridrow~=1
                error('Make sure that at least reac_id and reac_eq are located in the header line');
            end
        end
    end

    %%  ================== add everything to model ========================
    if cnap.has_gui
        % add map
        if ~isempty(map.map_path)
            try
                cd(cnap.path);
                cnap = CNAaddMap(cnap,map.map_path,map.map_name);
                cd(cnan.cnapath);
            catch ME
                % if adding the new map fails, put reactions on generated map
                cd(cnan.cnapath);
                warning(['map could not be added. Reactions are either placed on existing or on generated map'...
                    ,char(10),ME.identifier,': ',ME.message]);
                for reac = reactions
                    if reac.flxMapNo == 0
                        reactions([reactions.reac_id]==reac.reac).flxMapNo = -1;
                    end
                end
            end
        end
        % Adjust reaction box - flux map mapping for every reaction
        % Iterate all reactions to add and convert flxMapNo to integer value.
        % If flxMapNo is a character array, search for map name in cnap
        % If there are reactions to be positioned on the new map, update
        % their flxMapNo-field accordingly (set it to the index of the
        % new map.)
        % If everything fails, place reactions on generated map
        for reac = reactions
            index = find(strcmp(strtrim({reactions.reac}'),reac.reac));

            % 1. if a name is indicated, try to find corresponding map index
            reactions(index).flxMapNo = find(strcmp(cnap.maps(:,1),reac.flxMapNo));

            % 2. if a number is indicated, treat different cases
            %if isempty(reactions(index).flxMapNo)
                % 2.1. take integer value mentioned in xls-sheet
                reactions(index).flxMapNo = str2double(reac.flxMapNo);

                % 2.2. update map index to fit the index of the newly added map
                if reactions(index).flxMapNo == 0
                    reactions(index).flxMapNo = find(strcmp(cnap.maps(:,1),map.map_name));
                    % 2.3. if map name or index doesn't occur in cna project, place
                    % on generated map
                elseif isempty(reactions(index).flxMapNo)
                    reactions(index).flxMapNo = -1;
                elseif reactions(index).flxMapNo > size(cnap.maps(:,1),1) || reactions(index).flxMapNo < 0
                    reactions(index).flxMapNo = -1;
                end
            %end
        end
    end

    % add Species
    n = 1;  % index only needed for error output
    try
        for spec = species
            n = find(strcmp(strtrim([species.spec_id]'),spec.spec_id));
            cnap = CNAaddSpeciesMFN(cnap,spec);
        end
    catch ME
        % species already exists in model
        if length(n) ~= 1
            error('you defined multiple species with the same ID, make sure every ID appears only once.')            
        elseif any(strcmp(cellstr(cnap.specID),char(species(n).spec_id)))
            warning(['The species id "',char(species(n).spec_id),'" already occurs in the model.',char(10),ME.identifier,': ',ME.message]);
            % species cannot be added for other reasons
        else
            error(['The species id "',char(species(n).spec_id),'" could not be added to the model.',char(10),ME.identifier,': ',ME.message]);
        end
        % rethrow(ME);
    end

    % add Reactions
    n = 0;  % index only needed for error output
    try
        for reac = reactions
            n = find(strcmp(strtrim({reactions.reac}'),reac.reac));
            cnap = CNAaddReactionMFN(cnap,reac);
        end
    catch ME
        % reaction already exists in model
        if any(strcmp(cnap.reacID,char(reactions(n).reac)))
            error(['The reaction id "',char(reactions(n).reac),'" already occurs in the model.',char(10),ME.identifier,': ',ME.message]);
            % reaction cannot be added for other reasons
        else
            error(['The reaction id "',char(reactions(n).reac),'" could not be added to the model.',char(10),ME.identifier,': ',ME.message]);
        end

        % error(['Please make sure the reaction ids you defined, dont already occur in the model. No changes have been made to the project.',char(10),ME.identifier,': ',ME.message]);
    end
    % Remove all text boxes that don't correspond to reactions that are
    % contained in the project
    if cnap.has_gui
        deleteInvalidTextBoxes(cnap);
    % usually not necessary, but sometimes the generated map is not updated
    % correctly if this function isn't called
    %
    % With less than 50 reaction on the generated map, also show reaction equation
        cnap = myCNAgenerateMap(cnap,sum(cnap.reacBoxes(:,5) == -1) < 50);
    end
    
end