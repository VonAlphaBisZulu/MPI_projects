function cnap = CNAgenerateMap(cnap)
global cnan;
cnapBU = cnap; % backup
%% Generates a reaction map automatically (temporarily)
% cnap:         CNA project object
%
% Automatically generate a map for all reactions in a cnap project
% whose reaction boxes are either allocated to an invalid map index or
% to the generated map. Reactions are displayed in a vertical,
% list-like manner and reaction partners are shown on the two sides of
% an arrow. The reaction boxes are placed automatically.
% Changes by this function affect the cna project object and the gui,
% but are only temporarily. A dialog box will open up that proposes to
% store changes to the project permanently.
%
% Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
% Jun 15, 2017

if ~isfield(cnap,'path')
    error('The generated maps are stored in cnap.path, but cnap.path is undefined. If you''re coming from a new project, please save the network first.');
end
if ~cnap.has_gui
    if ~exist([cnap.path '/app_para.m'],'file')
        warning('CNA project doesn''t have a map yet. Using standard settings for map creation.');
        if ~isfield(cnap,'net_var_name')
            cnap.net_var_name = cnap.path; 
        end
        cnap.local.errval = 0;
        cnap.reacBoxes(:,5) = -1;
        cnap.nummaps = 0;
        generatedMaps = [];
        cnap.maps = cell(0,2);
        cnap.figs = double.empty(0,6);
        cnap.color1 = [0.7 0.7 0.7];
        cnap.color2 = [0.5 0.5 1  ];
        cnap.color3 = [1.  0.2 0.2];
        cnap.color4 = [0.2 1   0.2];
        cnap.textColor = [0 0 0];
        cnap.specBoxColor = [0.9 0.9 0.5];
        cnap.macroSynthColor = [0 0 1];
        cnap.reacBoxWidth = [];
        cnap.reacBoxHeight = [];
        cnap.specBoxWidth = [];
        cnap.specBoxHeight = [];
        cnap.reacFontSize = [];
        cnap.specFontSize = [];
    else
        disp('Please load the network with gui to append map.');
        return;
    end
else
    generatedMaps = [find(strcmp(cnap.maps(:,1),'generated')); find(strcmp(cnap.maps(:,1),'_generated'))];
end
dimensions.reacBoxWidth = 0.05;
dimensions.reacBoxHeight = 0.018;
dimensions.specBoxWidth = 0.04;
dimensions.specBoxHeight = 0.014;
dimensions.reacFontSize = 10;
dimensions.specFontSize = 10;


%% ========= identify reactions and their educts and products =========
% Produkte und Edukte der Reaktion ermitteln
% Liste mit zwei Spalten für Edukte und Produkte der Reaktionen
str=repmat({'','',''},0,1);
% get indices of generated maps and find reactions that are either already
% placed on them or allocated on an invalic map
reacsAlreadyOnMap = [];
if ~isempty(generatedMaps)
    reacsAlreadyOnMap = cnap.reacBoxes(:,5) == generatedMaps;
end
if isempty(reacsAlreadyOnMap)
    reacsAlreadyOnMap = zeros(size(cnap.reacBoxes,1),1);
end
reacsOnInvalidMap = cnap.reacBoxes(:,5) > cnap.nummaps | cnap.reacBoxes(:,5)==-1;
if isempty(reacsOnInvalidMap)
    reacsOnInvalidMap = zeros(size(cnap.reacBoxes,2),1);
end
reacsToBePlacedOnNewMap = find(or(reacsOnInvalidMap , reacsAlreadyOnMap));

% abort if list of reactions to be placed on generated map is empty
if isempty(reacsToBePlacedOnNewMap)
    cnap = cnapBU;
    if any(strcmp(cnap.maps(:,1),'generated'))
        disp('There are no reactions to be placed on a generated map. Generated Map was removed.\n');
        cnap = CNAremoveMap(cnap,'generated');
        return;
    end
    disp('There are no reactions to be placed on a generated map. No changes were applied.\n');
    return;
end

scrsize=get_screen_size();
r_per_column = 55;
r_columns = 4;

if scrsize(4)/40 >= length(reacsToBePlacedOnNewMap)-3 % maximum number of reactions that can be placed below each other
    % Find educts and products
    %
    %       {'name1'    'educt 1'    'product1'}
    % str = {'name2'    'educt 2'    'product2'}
    %       {'name3'    'educt 3'    'product3'}
    %       {'name4'    'educt 4'    'product4'}
    for reacIndex = reacsToBePlacedOnNewMap'
        educts='';
        zw=find(cnap.stoichMat(:,reacIndex)<0);
        
        if~isempty(zw)
            name = strtrim(strrep(cnap.reacID(reacIndex,:),'_','\_')); % make sure underscores are not "interpreted"
        end
        
        if~isempty(zw)
            educts=[num2str(-cnap.stoichMat(zw(1),reacIndex)),' ',strrep(deblank(cnap.specID(zw(1),:)),'_','\_')];
            for j=2:length(zw)
                educts=[educts,' + ',num2str(-cnap.stoichMat(zw(j),reacIndex)),' ',strrep(deblank(cnap.specID(zw(j),:)),'_','\_')];
            end
        end
        
        
        products='';
        zw=find(cnap.stoichMat(:,reacIndex)>0);
        if~isempty(zw)
            products=[products,num2str(cnap.stoichMat(zw(1),reacIndex)),' ',strrep(deblank(cnap.specID(zw(1),:)),'_','\_')];
            for j=2:length(zw)
                products=[products,' + ',num2str(cnap.stoichMat(zw(j),reacIndex)),' ',strrep(deblank(cnap.specID(zw(j),:)),'_','\_')];
            end
        end
        
        str=[str; {name,educts,products}];
    end
    
    %% ============== create the new fluxmap with arrows ==================
    
    % create reaction map figure (1200x40 , invisible)    
    newMap = figure('Position',[scrsize(3)/2-600 scrsize(4)/2+200 , 1200 , 40*size(str,1)],...
        'visible','off','Renderer','painters');
    a=axes('Position',[0 0 1 1]);
    axis off;
    
    % draw arrows (as many as there are reactions)
    % first arrow
    pLeft       = [0 0];
    pRight      = [1 0];
    
    for n = 1:size(str,1)-1
        pLeft   = [pLeft; 0 -n];
        pRight  = [pRight; 1 -n];
    end
    
    headLength  = 0.1;
    headWidth   = 0.2;
    
    for n = 1:size(str,1)
        headLeftX = [pLeft(n,1),  pLeft(n,1)+headLength,  pLeft(n,1)+headLength];
        headLeftY = [pLeft(n,2),  pLeft(n,2)-headWidth,  pLeft(n,2)+headWidth];
        
        headRightX = [pRight(n,1),  pRight(n,1)-headLength,  pRight(n,1)-headLength];
        headRightY = [pRight(n,2),  pRight(n,2)-headWidth,  pRight(n,2)+headWidth];
        
        arrows(n).shaft = line([pLeft(n,1) pRight(n,1)],  [pLeft(n,2) pRight(n,2)],'Color',[0 0 0],'LineWidth',2);
        arrows(n).leftHead = patch(headLeftX,headLeftY,[0 0 0],'Visible','off');
        arrows(n).rightHead = patch(headRightX,headRightY,[0 0 0]);
        if cnap.reacMin(reacsToBePlacedOnNewMap(n)) < 0
            set(arrows(n).leftHead,'Visible','on');
        end
    end
    
    % adjust image so that text boxes appear on the figure
    % 1. put textboxes
    set(a,'XLim',[0 1]) ;
    set(a,'YLim',[-n+0.5 0.5]) ;
    for n = 1:size(str,1)
        textl(n) = text(-0.02,pLeft(n,2), char(strcat(str(n,1),':  ',str(n,2))),'HorizontalAlignment','right');
        textr(n) = text(1.02,pLeft(n,2), char(str(n,3)),'HorizontalAlignment','left');
    end
    % 2. adjust position to fit margin left and right
    if isnumeric(get(textl,'Extent'))
        leftExt = get(textl,'Extent');
    elseif iscell(get(textl,'Extent'))
        leftExt = cell2mat(get(textl,'Extent'));
    end
    leftExt = min(leftExt(:,1));
    if isnumeric(get(textr,'Extent'))
        rightExt = get(textr,'Extent');
    elseif iscell(get(textr,'Extent'))
        rightExt = cell2mat(get(textr,'Extent'));
    end
    rigthExt = min(1+leftExt-rightExt(:,3)-0.02);
    set(a,'Position',[-leftExt 0 rigthExt 1]) ;
    % 3. left-align reaction names
    if isnumeric(get(textl,'Extent'))
        leftExt = get(textl,'Extent');
    elseif iscell(get(textl,'Extent'))
        leftExt = cell2mat(get(textl,'Extent'));
    end
    leftExt = min(leftExt(:,1));
    for n = textl
        oldExtents = get(n,'Extent');
        position = get(n,'Position');
        position(1) = position(1) + leftExt - oldExtents(1);
        set(n,'Position', position);
        %n.Position(1) = n.Position(1) + leftExt - n.Extent(1);
    end
    
    %set(newMap,'Visible','on');
    axis off;
    % print map as pcx to file
    set(newMap,'Units','inches');
    format = get(newMap,'Position');
    set(newMap,'PaperUnits','inches');
    if verLessThan('matlab', '8.4')
        set(newMap,'PaperPosition',1.25*format);
        saveas(newMap,[cnap.path,'/_generated.pcx'])
    else
        set(newMap,'PaperPosition',format);
        print(newMap,strcat(cnap.path,'/_generated'),'-dpcx256','-r93');
    end
    
    close(newMap);
    
    
    % Set the map parameter of all reaction boxes to be placed on generated
    % map to -1
    for reacIndex = reacsToBePlacedOnNewMap'
        cnap.reacBoxes(reacIndex,5) = -1;
    end
    
    % Delete formerly generated maps
    generatedMaps = find(strcmp(cnap.maps(:,1),'generated') | strcmp(cnap.maps(:,1),'_generated'));
    
    % Delete all existing reaction boxes that are about to be placed on the new
    % map. Then remove the old generated map(s). (Usually there should only be
    % none or one)
    for generMapIndex = generatedMaps'
        for reacIndex = reacsToBePlacedOnNewMap'
            try
                delete(cnap.reacBoxes(reacIndex,4));
            catch
                % warning disabled, enable for debugging
                % warning('Old reaction box not found. No element was deleted')
            end
        end
        cnap = CNAremoveMap(cnap,generMapIndex);
    end
    
    %% =========== Add new fluxmap and place reaction boxes ===============
    
    % load new Map _generated. (When saved, the generated map will be renamed
    % to "generated")
    cnap = CNAaddMap(cnap,strcat(cnap.path,'/_generated.pcx'),'generated');
    cnap.local.boolGeneratedMapWasUpdated = 1;
    
    % put the reaction boxes on the new map
    % Section was copied from reacnewmaskeval, line 110ff (Jun 14, 2017)
    % and adapted
    y = 0;  % y position of the new map
    for reacIndex = reacsToBePlacedOnNewMap'
        % reaction box position
        cnap.reacBoxes(reacIndex,5) = cnap.nummaps;
        cnap.reacBoxes(reacIndex,2) = 630;
        cnap.reacBoxes(reacIndex,3) = 32+y;
        y = y+40;
    end
else
    r_idx_placed = 0;
    r_per_column = 55;
    r_columns = 4;
    r_per_map = r_per_column * r_columns;
    rowHeight = 14;
    colWidth = 280;
    scrsize=get_screen_size();
    for i = reacsToBePlacedOnNewMap(:)'
        cnap.reacBoxes(i,2) = 200 + colWidth*floor(mod(r_idx_placed,r_per_map)/r_per_column); % xpos
        cnap.reacBoxes(i,3) = rowHeight  +  rowHeight*+mod(mod(r_idx_placed,r_per_map),r_per_column); % ypos
        cnap.reacBoxes(i,5) = -ceil((r_idx_placed+1)/r_per_map);% map
        r_idx_placed = r_idx_placed+1;
    end
    maps = flip(unique(cnap.reacBoxes((cnap.reacBoxes(:,5)<0),5)))';
    for i = maps
        disp(['Creating map ' num2str(-i) ' of ' num2str(max(abs(maps)))]);
        figWidth(-i)  = 100 + max(cnap.reacBoxes(cnap.reacBoxes(:,5)==i,2));
        figHeight(-i) = max(cnap.reacBoxes(cnap.reacBoxes(:,5)==i,3));
        newMap(-i) = figure('Units','pixels','Position',[scrsize(3)/2-600 scrsize(4)/2+200 , figWidth(-i), figHeight(-i)],...
            'visible','off','Renderer','painters','Color','white');
        reacs_on_map = find(cnap.reacBoxes(:,5)==i)';
        rows = unique(cnap.reacBoxes(reacs_on_map,3));
        for j = 1:(length(rows)/2)
            rctg(-i,j) = annotation(newMap(-i),'rectangle','Color','none',...
                'FaceColor',[0.9 0.93 0.97],'Units','pixels','Position',[0 rows(j*2-1)+1 figWidth(-i) rowHeight]);
        end
        for j = reacs_on_map
            str = cnap.reacID(j,:);
            if length(strtrim(str))>23
                str = [str(1:21) '..'];
            end        
            str = strrep(str,'_','\_');
            r_text(-i,j) = annotation(newMap(-i),'textbox','Units','pixels','EdgeColor','none','Color','black',...
                'String',str,'Position',[cnap.reacBoxes(j,2)-200 ...
                                           figHeight(-i)-(cnap.reacBoxes(j,3)-1.5*rowHeight)    0 0]);
        end
        set(newMap(-i),'Units','inches');
        format = get(newMap(-i),'Position');
        set(newMap(-i),'PaperUnits','inches');
        filename{-i} = strcat(cnap.path,'/generated_map_',num2str(-i,'%.2i'), '.pcx');
        if verLessThan('matlab', '8.4')
            set(newMap(-i),'PaperPosition',1.25*format);
            saveas(newMap(-i),filename{-i})
        else
            set(newMap(-i),'PaperPosition',format/300);
            write_pcx(newMap(-i),filename{-i});
            %print(newMap(-i),filename{-i},'-dpcx24b','-r93');
        end
    end
    for i = maps
        dimensions.reacBoxWidth = round(57/figWidth(-i),3);
        dimensions.reacBoxHeight = round(13.5/figHeight(-i),3);
        reacs_on_map = cnap.reacBoxes(:,5) == i;
        cnap.reacBoxes(reacs_on_map,5) = cnap.nummaps+1;
        cnap.unsaved_changes = 1;
        cnap.local.boolGeneratedMapWasUpdated = 1;
        cnap = CNAaddMap(cnap,filename{-i},['generated_map_' num2str(-i,'%.2i')],dimensions,1);
    end
end
%% =========== Place reaction boxes ===============

% put the reaction boxes on the new map
% Section was copied from reacnewmaskeval, line 110ff (Jun 14, 2017)
% and adapted
disp('Placing reaction boxes');
for reacIndex = reacsToBePlacedOnNewMap'
    % declaration of handle
    cnan.open_projects.(cnap.net_var_name).gui.handles= struct;
    cnap = update_after_change(cnap);
    currmap = cnap.reacBoxes(reacIndex,5);
    % is visible/editable
    zw=cnap.reacBoxes(reacIndex,6);
    fig = cnap.figs(currmap,:);
    % make box
    zw1=uicontrol('Style', 'edit','Parent',fig(1), 'String', '###',...
        'Units','normalized','HorizontalAlignment','left','BackgroundColor',cnap.color1,'ForegroundColor',cnap.textColor,'TooltipString',cnap.reacID(reacIndex,:));
    set(zw1, 'ButtonDownFcn', {@execute_callback, cnap.net_var_name,...
        {'check_for_right_click', 'reaceditmask'}, {'reacenr', reacIndex}});
    % save handle
    cnap.reacBoxes(reacIndex,4)=zw1;
    % adjust "zoom"
    place_box(fig,zw1,...
        cnap.reacBoxes(reacIndex,2),...
        cnap.reacBoxes(reacIndex,3),...
        cnap.reacFontSize(currmap),...
        cnap.reacBoxWidth(currmap),...
        cnap.reacBoxHeight(currmap));
    % put default rate in textbox
    if isnan(cnap.reacDefault(reacIndex))
        set(zw1,'String','#');
    else
        set(zw1,'String',num2str(cnap.reacDefault(reacIndex)));
    end
    % is visible/editable
    if(zw==2)
        set(zw1,'Style', 'text');
    elseif(zw==3)  %%non-visible
        set(zw1,'Visible','off');
    end
end
% delete reaction boxes that are not acessible through the handles in cnap.reacBoxes(:,4)
deleteInvalidTextBoxes(cnap);
set(cnap.figs(:,1),'Visible','on');
saveMapCfgDiag(cnap);

function write_pcx(fig,filename)
    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256); % <= 65536
    imwrite(A,map,filename);
end

function place_box(fig,handle,xp,yp,fontsize,box_width,box_height) % copied from zoom_single_box
    zzz=get(fig(1),'CurrentAxes');
    if(~fig(4))
    xxx=get(zzz,'XLIM');
    yyy=get(zzz,'YLIM');
    xxxl=xxx(2)-xxx(1);
    yyyl=yyy(2)-yyy(1);
    xp=(xp-xxx(1))/xxxl;
    yp=(yyy(2)-yp)/yyyl;
    box_width=box_width*fig(3)/xxxl;
    box_height=box_height*fig(2)/yyyl;
    fontsize=fontsize*fig(3)/xxxl;
    end
    pos=[xp yp box_width box_height];
    set(handle,'Position',pos,'FontSize',fontsize);
    end
end