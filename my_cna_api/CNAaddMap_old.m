function cnap = CNAaddMap( cnap, imagepath, mapname, dimensions, nodisp)
    global cnan;
    %% Adds a new map to the project (temporarily)
    %
    % Usage:  cnap = CNAaddMap( cnap, 'relative/path/to/file.pcx', 'mapname' (, dimensions) )
    %
    % Input:
    %
    %   cnap:         CNA project object
    %   imagepath:    [char]   relative path to image (from cnapath / CNA root)
    %   mapname:      [char]   name of the map to be saved at cnap.maps (also in header of
    %                        the map window)
    %   dimensions:   (struct) structure with fields that describe layout parameters:
    %                           (double) reacBoxWidth  (d) reacBoxHeight  (d) specBoxWidth
    %                           (d)      specBoxHeight (d) reacFontSize   (d) specFontSize
    % 
    % Output: 
    %   
    %   cnap:         CNA project object
    %
    % Adds a new map to the project (temporarily). Permanent modifications to
    % the poject can either be made by saving the map configuration through the
    % dialog or by calling CNAsaveProjectWithNewMap
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 13, 2017

    %% ======================= process input ==============================
    % cut unnecessary ./ from path
    if all(imagepath(1,2)=='./')
        imagepath = imagepath(3:end);
    end
    % check if map image file exists and has the correct type (png or pcx)
    if exist(imagepath,'file') == 0
        error('Path to image invalid');
    else
        [filepath,filename,ext] = fileparts(imagepath);
        if strcmp(ext,'.pcx') && strcmp(ext,'.png')
            error('Bad image type. Choose pcx or png');
        end
    end

    % if no name is chosen, map is named after file
    if nargin == 2
        mapname = filename;
    end
    if nargin < 4
        dimensions = [];
    end
    if nargin <5
        nodisp = 0;
    end
    
    % if image file is in project folder no path needs to be indicated;
    % only the filename
    if strcmp(filepath,cnap.path)
        relimgpath = [filename ext];
        % else indication to go one level up is needed, because the path to
        % image is always referred to from the cna project folder.
    else
        relimgpath = [filepath '/' filename ext];
    end

    % check if map name is already defined in cna project
    if ~isempty(find(strcmp(cnap.maps(:,1),mapname),1))
        error ('Map name already in project');
    end
    % get image dimensions from file
    img.map    = imread(imagepath);
    img.height = size(img.map,1);
    img.width  = size(img.map,2);

    %% ================== modify project object ===========================
    % cna project object parameters are changed to contain new map

    % set cnap parameters related to the flux maps, to add file
    cnap.maps = [cnap.maps ; {char(mapname) char(relimgpath)}];
    if ~isempty(cnap.figs)
        cnap.figs = [cnap.figs ; [0 img.height img.width cnap.figs(1,4:end)]];
    else
        cnap.figs = [0 img.height img.width 0 0 0 0 0 0];
    end
    cnap.nummaps = cnap.nummaps+1;

    % set style parameters for the map
    % if no dimensions for the new map are given, reaction box sizes are
    % calculated automatically or derived from existing maps.
    % Those styke parameters are stored in vectors that are propagated
    if isempty(dimensions)
        cnap.reacBoxWidth  = [cnap.reacBoxWidth 55/img.width];
        cnap.reacBoxHeight = [cnap.reacBoxHeight 25/img.height];
        cnap.specBoxWidth  = [cnap.specBoxWidth 55/img.width];
        cnap.specBoxHeight = [cnap.specBoxHeight 25/img.height];
        cnap.reacFontSize  = [cnap.reacFontSize cnap.reacFontSize(1)]; % Inherited from first map
        cnap.specFontSize  = [cnap.specFontSize cnap.specFontSize(1)];
        % if dimensions structure is passed, it is verified whether all fields
        % needed are present
    elseif isfield(dimensions,'reacBoxWidth') && ...
            isfield(dimensions,'reacBoxHeight') && ...
            isfield(dimensions,'specBoxWidth') && ...
            isfield(dimensions,'specBoxHeight') && ...
            isfield(dimensions,'reacFontSize') && ...
            isfield(dimensions,'specFontSize')
        cnap.reacBoxWidth  = [cnap.reacBoxWidth dimensions.reacBoxWidth];
        cnap.reacBoxHeight = [cnap.reacBoxHeight dimensions.reacBoxHeight];
        cnap.specBoxWidth  = [cnap.specBoxWidth dimensions.specBoxWidth];
        cnap.specBoxHeight = [cnap.specBoxHeight dimensions.specBoxHeight];
        cnap.reacFontSize  = [cnap.reacFontSize dimensions.reacFontSize];
        cnap.specFontSize  = [cnap.specFontSize dimensions.specFontSize];
    else
        error('When you call this function and pass the dimension argument, make sure it is a sructure that contains all required fields');
    end
    
    %% =================== load fluxmap in window =========================

    % declaration: project is set as open project in cnan
%     cnan.open_projects.(cnap.net_var_name).gui.handles= struct;
%     cd(cnan.cnapath);
%     cd(cnap.path);
    try
        IndexLastMapAdded=max(cnap.nummaps);
%         % set temporary map-index-indicator in cna project
%         cnap.local.mapnr=IndexLastMapAdded;
% 
%         % load image
%         % visibility will be turned on at the end
%         cnap.figs(IndexLastMapAdded, 1)= figure('Visible', 'off');
%         cnap.figs(IndexLastMapAdded, 4)= 0;   % Typ:Matlab-Figure
%         [imagedata,colmap] = imread(imagepath);
%         colormap(colmap);
%         image(imagedata);
%         zw=size(imagedata);
%         cnap.figs(IndexLastMapAdded,2:3)=zw(1:2);
%         set(get(cnap.figs(IndexLastMapAdded,1),'Children'),'Position',[0,0,1,1]);
%         set(get(cnap.figs(IndexLastMapAdded,1),'Children'),'Visible', 'off');
% 
%         % set up toolbar
%         cnap = install_toolbar(cnap);
%         set(cnap.figs(IndexLastMapAdded,1),'NumberTitle','off','Name',['CellNetAnalyzer: ',cnap.maps{IndexLastMapAdded,1}, ' (', cnap.net_var_name, ')'],'MenuBar','none');
%         set(cnap.figs(IndexLastMapAdded,1),'Units','pixels');
%         set(cnap.figs(IndexLastMapAdded,1),'CloseRequestFcn', ['closeDialog2(''UserData'',''', cnap.net_var_name, ''');']);
        cnan.open_projects.(cnap.net_var_name).gui.handles= struct;
        cd(cnan.cnapath);
        cd(cnap.path);
        cnap = flux_init(cnap, ~cnap.has_gui , IndexLastMapAdded);

        %% ==== adapt position and resolution ====
        % complete copy from optisize function (Jun 13, 2017)

        set(cnap.figs(IndexLastMapAdded,1),'Units','pixels');
        scrsize=get_screen_size();
        zw1=cnap.figs(IndexLastMapAdded,3);
        zw2=cnap.figs(IndexLastMapAdded,2);
        if(min(scrsize(3)-zw1-320,scrsize(4)-zw2-90)>0)
            set(cnap.figs(IndexLastMapAdded,1),'Position',[280 50 zw1 zw2]);
        elseif(min(scrsize(3)-zw1,scrsize(4)-zw2)>0)
            set(cnap.figs(IndexLastMapAdded,1),'Position',[1 1 zw1 zw2]);
        else
            width= scrsize(3) - 80;
            height= scrsize(4) - 80;
            f_x= width/zw1;
            f_y= height/zw2;
            if f_x < f_y
                f= f_x;
            else
                f= f_y;
            end
            set(cnap.figs(IndexLastMapAdded,1), 'Position', [10 10 zw1*f zw2*f]);
            fprintf('Hint: Network map %d can not be displayed in full resolution, activating zoom tools.\n\n', 1i);
            figure(cnap.figs(IndexLastMapAdded, 1)); %A# make sure the gcf() calls in the following two functions work correctly
            cnap= clip_axis_to_window(cnap);
            cnap= zoomonoff(cnap);
        end
        pos= get(cnap.figs(IndexLastMapAdded,1), 'Position'); %A# get actual position
        user_data.width= pos(3);
        user_data.height= pos(4);
        set(cnap.figs(IndexLastMapAdded,1), 'UserData', user_data);

        % make visible
        if ~nodisp
            figure(cnap.figs(1,1));
            set(cnap.figs(:, 1), 'Visible', 'on');
        end
    catch EM
        warning(['New Map could not be loaded.',char(10),EM.identifier,': ',EM.message]);
        cnap.local.errval= 1;
    end

    cd(cnan.cnapath);
    if cnap.local.errval == 0
        assignin('base', cnap.net_var_name, cnap);
        fprintf('New map created: %s\n', mapname);
        % in case of error remove project from cnan
    else
%         cnan.open_projects= rmfield(cnan.open_projects, cnap.net_var_name);
    end
    %% ====== dialog to save changes in map configuration permanently =====
    if ~nodisp
        saveMapCfgDiag(cnap);
    end
end