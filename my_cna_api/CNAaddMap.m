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
    % Copied from flux_init: line 76ff  (Jun 13, 2017) and only
    % slightly adapted
    %
    % From here on functions are called to load the map in a window, set up
    % figure handles, load the menu and menu bar, etc.

    % declaration: project is set as open project in cnan
    cnan.open_projects.(cnap.net_var_name).gui.handles= struct;
    cd(cnan.cnapath);
    cd(cnap.path);
    try
        % only for metabolic networks
        if cnap.type == 1
            %% ==== load image and toolbar ====
            IndexLastMapAdded=max(cnap.nummaps);
            % set temporary map-index-indicator in cna project
            cnap.local.mapnr=IndexLastMapAdded;

            % load image
            % visibility will be turned on at the end
            cnap.figs(IndexLastMapAdded, 1)= figure('Visible', 'off');
            cnap.figs(IndexLastMapAdded, 4)= 0;   % Typ:Matlab-Figure
            [imagedata,colmap] = imread(cnap.maps{IndexLastMapAdded,2});
            colormap(colmap);
            image(imagedata);
            zw=size(imagedata);
            cnap.figs(IndexLastMapAdded,2:3)=zw(1:2);

            % set up toolbar
            cnap = install_toolbar(cnap);
            set(cnap.figs(IndexLastMapAdded,1),'NumberTitle','off','Name',['CellNetAnalyzer: ',cnap.maps{IndexLastMapAdded,1}, ' (', cnap.net_var_name, ')'],'MenuBar','none');
            set(cnap.figs(IndexLastMapAdded,1),'Units','pixels');
            set(cnap.figs(IndexLastMapAdded,1),'CloseRequestFcn', ['closeDialog2(''UserData'',''', cnap.net_var_name, ''');']);
            set(gca,'Position',[0,0,1,1]);
            set(gca,'Visible','off');

            %% ==== load reaction boxes to map (they should not yet exist) ====
            % Still reactions can already be assigned to a high map index and
            % now appear on the newly added map

            % GUIS for Reactions (should usually not exist at that point)
            if(cnap.numr>0)
                zw=find(cnap.reacBoxes(:,5)==cnap.local.mapnr);

                for li=1:length(zw)
                    if(cnap.figs(IndexLastMapAdded,4))
                        xp=cnap.reacBoxes(zw(li),2);
                        yp=cnap.reacBoxes(zw(li),3);
                    else
                        xp=cnap.reacBoxes(zw(li),2)/cnap.figs(cnap.local.mapnr,3);
                        yp=(cnap.figs(cnap.local.mapnr,2)-cnap.reacBoxes(zw(li),3))/cnap.figs(cnap.local.mapnr,2);
                    end

                    if(cnap.reacBoxes(zw(li),6)==1)
                        cnap.reacBoxes(zw(li),4)=uicontrol('Style', 'edit', 'String', '#','Units','normalized','Position',[xp,yp,cnap.reacBoxWidth(cnap.local.mapnr),cnap.reacBoxHeight(cnap.local.mapnr)],'FontSize',cnap.reacFontSize(cnap.local.mapnr),'HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',deblank(cnap.reacID(zw(li),:)));

                    else
                        cnap.reacBoxes(zw(li),4)=uicontrol('Style','text', 'String', '#','Units','normalized','Position',[xp,yp,cnap.reacBoxWidth(cnap.local.mapnr),cnap.reacBoxHeight(cnap.local.mapnr)],'FontSize',cnap.reacFontSize(cnap.local.mapnr),'HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',deblank(cnap.reacID(zw(li),:)));
                    end
                    set(cnap.reacBoxes(zw(li),4), 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
                        {'check_for_right_click', 'reaceditmask'}, {'reacenr', zw(li)}});
                end
            else
                cnap.reacBoxes= zeros(0, 6);
            end


            % GUIS for Macromolecules (should usually not exist at that point)
            if(cnap.nummac>0)
                zw=find(cnap.macroBoxes(:,5)==cnap.local.mapnr);

                for li=1:length(zw)
                    if(cnap.figs(cnap.local.mapnr,4))
                        xp=cnap.macroBoxes(zw(li),2);
                        yp=cnap.macroBoxes(zw(li),3);
                    else
                        xp=cnap.macroBoxes(zw(li),2)/cnap.figs(cnap.local.mapnr,3);
                        yp=(cnap.figs(cnap.local.mapnr,2)-cnap.macroBoxes(zw(li),3))/cnap.figs(cnap.local.mapnr,2);
                    end

                    if(cnap.macroBoxes(zw(li),6)==1)
                        cnap.macroBoxes(zw(li),4)=uicontrol('Style', 'edit', 'String', '#','Units','normalized','Position',[xp,yp,cnap.specBoxWidth(cnap.local.mapnr),cnap.specBoxHeight(cnap.local.mapnr)],'BackgroundColor',cnap.specBoxColor,'ForegroundColor',cnap.textColor,'FontSize',cnap.specFontSize(cnap.local.mapnr),'HorizontalAlignment','Left','TooltipString',deblank(cnap.macroLongName(zw(li),:)));
                    else
                        cnap.macroBoxes(zw(li),4)=uicontrol('Style', 'text', 'String', '#','Units','normalized','Position',[xp,yp,cnap.specBoxWidth(cnap.local.mapnr),cnap.specBoxHeight(cnap.local.mapnr)],'BackgroundColor',cnap.specBoxColor,'ForegroundColor',cnap.textColor,'FontSize',cnap.specFontSize(cnap.local.mapnr),'HorizontalAlignment','Left','TooltipString',deblank(cnap.macroLongName(zw(li),:)));
                    end
                    set(cnap.macroBoxes(zw(li),4), 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
                        {'check_for_right_click', 'maceditmask'}, {'macenr', zw(li)}});
                end
            else
                cnap.macroBoxes= zeros(0, 6);
            end

            % GUIS for Macromolecule Synthesis (should usually not exist at that point)
            if(cnap.nummacsynth>0)
                zw=find(cnap.macroSynthBoxes(:,6)==cnap.local.mapnr);

                for li=1:length(zw)
                    if(cnap.figs(cnap.local.mapnr,4))
                        xp=cnap.macroSynthBoxes(zw(li),3);
                        yp=cnap.macroSynthBoxes(zw(li),4);
                    else
                        xp=cnap.macroSynthBoxes(zw(li),3)/cnap.figs(cnap.local.mapnr,3);
                        yp=(cnap.figs(cnap.local.mapnr,2)-cnap.macroSynthBoxes(zw(li),4))/cnap.figs(cnap.local.mapnr,2);
                    end

                    cnap.macroSynthBoxes(zw(li),5)=text('String', '#','Units','normalized','Position',[xp,yp],'Color',cnap.macroSynthColor,'FontSize',cnap.reacFontSize(cnap.local.mapnr));
                end
            else
                cnap.macroSynthBoxes= zeros(0, 6);
            end

            %% ==== add menu to window ====
            i = IndexLastMapAdded;
            % exact copy from inter_init, line 178ff (Jun 13, 2017)
            
            %% Menue
            cnap.figmenu(i,1)=uimenu(cnap.figs(i,1),'Label','Network');
            uimenu(cnap.figmenu(i,1),'Label','Network composer ...','Callback',[cnap.net_var_name, '=netcomposer(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,1),'Label','Save network ...','Callback',[cnap.net_var_name, '=save_net(', cnap.net_var_name, ');'],'Separator','off');
            exfigm=uimenu(cnap.figmenu(i,1),'Label','Export ...','Separator','off');
            uimenu(exfigm,'Label','Export stoichiometric matrix ...','Callback',[cnap.net_var_name, '=save_stoichmat_menu(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(exfigm,'Label','Export in SBML format ...','Callback',[cnap.net_var_name, '=export2sbml(', cnap.net_var_name, ');']);
            uimenu(exfigm,'Label','Export in COBRA toolbox format ...','Callback',[cnap.net_var_name, '=export2cobra(', cnap.net_var_name, ');']);
            uimenu(exfigm,'Label','Export in METATOOL format ...','Callback',[cnap.net_var_name, '=export2metatool(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,1),'Label','Import SBML model ...','Callback',['sbml_cna_import(', cnap.net_var_name, ');'],'Separator','off');

            uimenu(cnap.figmenu(i,1),'Label','Element selector ...','Callback',[cnap.net_var_name, '=open_element_selector(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,1),'Label','Find reaction ...','Callback',{@findReacDiag,cnap},'Separator','off');
            uimenu(cnap.figmenu(i,1),'Label','Show all reaction equations','Callback',[cnap.net_var_name, '=show_reactions(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(cnap.figmenu(i,1),'Label','Show names of network elements','Callback',[cnap.net_var_name, '=show_element_names(', cnap.net_var_name, ');']);
            %	uimenu(cnap.figmenu(i,1),'Label','Show constraints','Callback',[cnap.net_var_name, '=show_constraints(', cnap.net_var_name, ');'],'Separator','off');
            %	uimenu(cnap.figmenu(i,1),'Label','Set constraints','Callback',[cnap.net_var_name, '=get_lp_paras(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,1),'Label','Set epsilon and flux display format ...','Callback',[cnap.net_var_name, '=set_epsilon_flux_format(', cnap.net_var_name, ');'],'Separator','on');

            cnap.figmenu(i,2)=uimenu(cnap.figs(i,1),'Label','Map');
            uimenu(cnap.figmenu(i,2),'Label','Generate standard maps ...','Callback',['start_generate_maps(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,2),'Label','Set original map size','Callback',[cnap.net_var_name, '=optisize(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,2),'Label','Zoom tools on/off','Callback',[cnap.net_var_name, '=zoomonoff(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,2),'Label','Save map in graphics file ...','Callback',['save_map(',num2str(cnap.figs(i,1)),');']);

            cnap.figmenu(i,3)=uimenu(cnap.figs(i,1),'Label','Scenario');
            uimenu(cnap.figmenu(i,3),'Label','Clear all values','Callback',[cnap.net_var_name, '=clear_all(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,3),'Label','Reset last scenario','Callback',[cnap.net_var_name, '=set_last_scenario(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,3),'Label','Set default scenario','Callback',[cnap.net_var_name, '=set_default_values(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,3),'Label','Save scenario (fluxes only) ...','Callback',[cnap.net_var_name, '=save_values(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,3),'Label','Save scenario (fluxes, flux bounds, objective) ...','Callback',[cnap.net_var_name, '=save_values_all(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(cnap.figmenu(i,3),'Label','Load scenario ...','Separator','off','Callback',...
                {@execute_callback, cnap.net_var_name, {'get_file_name', 'get_rates', 'load_values', 'show_known_rates', 'highlight_values'},...
                {'dialog_title', 'Load file', 'filter_spec', '*.val'}});
            uimenu(cnap.figmenu(i,3),'Label','Highlight values (on/off)','Callback',[cnap.net_var_name, '=highlight_values(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,3),'Label','Highlight values (heatmap)','Callback',[cnap.net_var_name, '=heatmap_values(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,3),'Label','Show flux values in bar chart','Callback',[cnap.net_var_name, '=fluxbar(', cnap.net_var_name, ');']);

            cnap.figmenu(i,4)=uimenu(cnap.figs(i,1),'Label','Clipboard');
            uimenu(cnap.figmenu(i,4),'Label','Copy values to flux clipboard','Callback',[cnap.net_var_name, '=flux2clip(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(cnap.figmenu(i,4),'Label','Paste values from flux clipboard','Callback',[cnap.net_var_name, '=clip2map(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,4),'Label','Arithmetic operations ...','Callback',[cnap.net_var_name, '=arithmetic_ops_menu(', cnap.net_var_name, ');']);

            cnap.figmenu(i,5)=uimenu(cnap.figs(i,1),'Label','Analysis');
            %uimenu(cnap.figmenu(i,4),'Label','Show stoichiometric matrix','Callback',[cnap.net_var_name, '=plotmatstart(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,5),'Label','Basic network properties','Callback',[cnap.net_var_name, '=network_prop(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(cnap.figmenu(i,5),'Label','Conservation relations ...','Callback',[cnap.net_var_name, '=startconsrelations(', cnap.net_var_name, ');']);
            graphmen=uimenu(cnap.figmenu(i,5),'Label','Graph properties');
            uimenu(graphmen,'Label','Connectivity histogram','Callback',[cnap.net_var_name, '=metconhist(', cnap.net_var_name, ');']);
            uimenu(graphmen,'Label','Graph-theoretical path lengths ...','Callback',[cnap.net_var_name, '=startgraphpath(', cnap.net_var_name, ');'],'Separator','off');
            pathmen=uimenu(cnap.figmenu(i,5),'Label','Elementary flux modes/vectors','Separator','on');
            uimenu(pathmen,'Label','Compute ...','Callback',[cnap.net_var_name, '=startelmodes(', cnap.net_var_name, ');']);
            uimenu(pathmen,'Label','Load ...','Callback',[cnap.net_var_name, '=load_elmodes(', cnap.net_var_name, ');']);
            uimenu(pathmen,'Label','Show ...','Callback',[cnap.net_var_name, '=show_elmodes(', cnap.net_var_name, ');']);
            designmen=uimenu(cnap.figmenu(i,5),'Label','Network/strain design');
            mcsmen=uimenu(designmen,'Label','Minimal cut sets');
            uimenu(mcsmen,'Label','Compute ...','Callback',[cnap.net_var_name, '=start_cutsets_calc(', cnap.net_var_name, ');']);
            uimenu(mcsmen,'Label','Load ...','Callback',[cnap.net_var_name, '=load_cutsets(', cnap.net_var_name, ');']);
            uimenu(mcsmen,'Label','Show ...','Callback',[cnap.net_var_name, '=show_cutsets(', cnap.net_var_name, ');']);
            uimenu(designmen,'Label','CASOP ...', 'Callback', [cnap.net_var_name, '=start_CASOP(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,5),'Label','Check feasibility of flux scenario','Callback',[cnap.net_var_name, '=feasibility(', cnap.net_var_name, ');'],'Separator','on');
            mfamen=uimenu(cnap.figmenu(i,5),'Label','Metabolic flux analysis');
            uimenu(mfamen,'Label','Classify rates (balanceability/calculability) ','Callback',[cnap.net_var_name, '=classify(', cnap.net_var_name, ');']);
            uimenu(mfamen,'Label','Flux analysis ...','Callback',[cnap.net_var_name, '=determine_rates(', cnap.net_var_name, ');']);
            uimenu(mfamen,'Label','Sensitivity analysis ...','Callback',[cnap.net_var_name, '=sens_start(', cnap.net_var_name, ');']);
            uimenu(cnap.figmenu(i,5),'Label','Flux variability analysis','Callback',[cnap.net_var_name, '=flux_variability(', cnap.net_var_name, ');'],'Separator','off');
            uimenu(cnap.figmenu(i,5),'Label','Phase plane analysis ...','Callback',[cnap.net_var_name, '=phaseplane_plot(', cnap.net_var_name, ');'],'Separator','off');
            yieldmen=uimenu(cnap.figmenu(i,5),'Label','Yield analysis');
            uimenu(yieldmen,'Label','Maximize yield ...','Callback',[cnap.net_var_name, '=yield_optimization(', cnap.net_var_name, ');']);
            uimenu(yieldmen,'Label','2D yield space plot ...','Callback',[cnap.net_var_name, '=yieldspace_plot(', cnap.net_var_name, ');']);
            fbamen=uimenu(cnap.figmenu(i,5),'Label','Flux balance analysis (FBA)');
            uimenu(fbamen,'Label','Show objective function','Callback',[cnap.net_var_name, '=show_obj_func(', cnap.net_var_name, ');']);
            uimenu(fbamen,'Label','Flux optimization (FBA)','Callback',[cnap.net_var_name, '=optimize_flux(', cnap.net_var_name, ',0);']);
            uimenu(fbamen,'Label','Parsimonious FBA','Callback',[cnap.net_var_name, '=optimize_flux(', cnap.net_var_name, ',1);']);
            uimenu(cnap.figmenu(i,5),'Label','Net conversion of external metabolites','Callback',[cnap.net_var_name, '=net_conversion(', cnap.net_var_name, ');'],'Separator','on');
            uimenu(cnap.figmenu(i,5),'Label','In/Out fluxes at a metabolite ...','Callback',[cnap.net_var_name, '=in_out_fluxes(', cnap.net_var_name, ');']);

            cnap.figmenu(i,6)=uimenu(cnap.figs(i,1),'Label','Info');
            uimenu(cnap.figmenu(i,6),'Label','CellNetAnalyzer info ...','Callback','show_info');

            set(cnap.figs(IndexLastMapAdded,1),'KeyPressFcn',{@findReacDiag , cnap});
            set(cnap.reacBoxes(cnap.reacBoxes(:,5)==IndexLastMapAdded,4),'KeyPressFcn',{@findReacDiag , cnap});

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
        cnan.open_projects= rmfield(cnan.open_projects, cnap.net_var_name);
    end
    %% ====== dialog to save changes in map configuration permanently =====
    if ~nodisp
        saveMapCfgDiag(cnap);
    end
end