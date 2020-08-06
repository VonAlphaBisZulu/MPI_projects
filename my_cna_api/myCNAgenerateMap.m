function [cnap,err] = myCNAgenerateMap(cnap,show_eq)
%
% CellNetAnalyzer API function 'CNAgenerateMap'
% ---------------------------------------------
% CNA API function that creates maps for a given CNA mass-flow project listing 
% all reaction names. These maps can be used for a GUI-based analysis of the project.
%
% Usage: [cnap,err] = CNAgenerateMap(cnap,show_eq)
%
% On the generated maps, reaction names (optionally with reaction equations) 
% are displayed  in a vertical, list-like manner. The reaction boxes are 
% placed automatically. Several maps will be created depending on the number 
% of reactions and the maps (generated_map_01.pcx; generated_map_02.pcx; ....)
% are stored in the project folder.
%
% Note: elements of type "biomass constituents" will not be considered. 
%
% Inputs: 
% 	cnap:  CNA mass-flow project 
%
%       show_eq (optional; default: 0): whether the reaction equations are to
%                     be shown behind each reaction name.     
%
% Output:
%
% 	cnap: the new network project with each reaction text box placed on 
%             its position on the generate maps 
%
% 	err:  indicates, whether the an error occured during generation of the maps
%	      (err=1) or not (err=0).  

% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2019 by Steffen Klamt and Axel von Kamp and Philipp Schneider.
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


global cnan;
cnapBU = cnap; % backup
err=1;
if nargin < 2
    show_eq = 0;
end

str=repmat({'','',''},0,1);
reacs_on_generated_map = cell2mat(CNAgetGenericReactionData_as_array(cnap,'on_generated_map'));
if any(reacs_on_generated_map)
    cnap = CNAremoveMap(cnap,unique(cnap.reacBoxes(reacs_on_generated_map == 1,5)),1);
end
cnap.reacBoxes(reacs_on_generated_map == 1,5) = -1;
reacsToBePlacedOnNewMap = find(cnap.reacBoxes(:,5) == -1);
cnap = CNAsetGenericReactionData_with_array(cnap,'on_generated_map',num2cell(double(cnap.reacBoxes(:,5) == -1)));

scrsize=get_screen_size();

r_idx_placed = 0;
r_per_column = 55;
rowHeight = 14;

if show_eq
    colWidth = 560;
    % depending on the screen resolution, create up to 3 columns
    r_columns = min([3,floor(scrsize(3)/colWidth*0.75),ceil(length(reacsToBePlacedOnNewMap)/r_per_column)]);
    r_per_map = r_per_column * r_columns;
    for i = reacsToBePlacedOnNewMap(:)'
        mapnr=cnap.nummaps+ceil((r_idx_placed+1)/r_per_map);
        cnap.reacBoxes(i,2) = colWidth-80 + colWidth*floor(mod(r_idx_placed,r_per_map)/r_per_column); % xpos
        cnap.reacBoxes(i,3) = rowHeight  +  rowHeight*+mod(mod(r_idx_placed,r_per_map),r_per_column); % ypos
        cnap.reacBoxes(i,5) = mapnr;% map
        r_idx_placed = r_idx_placed+1;
        str{i} = getEQ(cnap,i);
        if length(strtrim(str{i}))>(colWidth/9.3)
            str{i} = [str{i}(1:round(colWidth/9.3)-2) '..'];
        end
        str{i} = strrep(str{i},'_','\_');
    end
else
    % depending on the screen resolution, create up to 4 columns
    colWidth = 280;
    r_columns = min(4,floor(scrsize(3)/colWidth*0.75));
    r_per_map = r_per_column * r_columns;
    for i = reacsToBePlacedOnNewMap(:)'
        mapnr=cnap.nummaps+ceil((r_idx_placed+1)/r_per_map);
        cnap.reacBoxes(i,2) = colWidth-80 + colWidth*floor(mod(r_idx_placed,r_per_map)/r_per_column); % xpos
        cnap.reacBoxes(i,3) = rowHeight  +  rowHeight*+mod(mod(r_idx_placed,r_per_map),r_per_column); % ypos
        cnap.reacBoxes(i,5) = mapnr;% map
        r_idx_placed = r_idx_placed+1;
        str{i} = cnap.reacID(i,:);
        if length(strtrim(str{i}))>23
            str{i} = [str{i}(1:21) '..'];
        end
        str{i} = strrep(str{i},'_','\_');
    end
end
newmaps = (cnap.nummaps+1):mapnr;
for i = newmaps
    disp(['Creating map ' num2str(i-cnap.nummaps) ' of ' num2str(mapnr-cnap.nummaps)]);
    figWidth(i)  = 100 + max(cnap.reacBoxes(cnap.reacBoxes(:,5)==i,2));
    figHeight(i) = max(cnap.reacBoxes(cnap.reacBoxes(:,5)==i,3));
    newMap = figure('Units','pixels','Position',[0, 0 , figWidth(i), figHeight(i)],...
        'Visible','off','Renderer','painters','Color','white');
    reacs_on_map = find(cnap.reacBoxes(:,5)==i)';
    rows = unique(cnap.reacBoxes(reacs_on_map,3));
    for j = 1:(length(rows)/2)
        rctg(i,j) = annotation(newMap,'rectangle','Color','none',...
            'FaceColor',[0.9 0.93 0.97],'Units','pixels','Position',[0 rows(j*2-1)+1 figWidth(i) rowHeight]);
    end
    for j = reacs_on_map
        r_text(i,j) = annotation(newMap,'textbox','Units','pixels','EdgeColor','none','Color','black',...
            'String',str{j},'Position',[cnap.reacBoxes(j,2)-colWidth+80 ...
            figHeight(i)-(cnap.reacBoxes(j,3)-1.5*rowHeight) colWidth 0]);
    end
%     cnap.figs(i,1:3) = [ newMap , figHeight(i), figWidth(i) ];
    filename{i} = strcat(cnap.path,'/generated_map_',num2str(i-cnap.nummaps,'%.2i'), '.png');
    frame = getframe(newMap);
    imwrite(frame.cdata,filename{i},'PNG');
    
    dimensions(i).reacBoxWidth = round(1000*57/figWidth(i))/1000;
    dimensions(i).reacBoxHeight = round(1000*13.5/figHeight(i))/1000;
    dimensions(i).specBoxWidth = round(1000*57/figWidth(i))/1000;
    dimensions(i).specBoxHeight = round(1000*13.5/figHeight(i))/1000;
    dimensions(i).reacFontSize  = 10; % Inherited from first map
    dimensions(i).specFontSize  = 10;
end
for i = newmaps
    reacs_on_map = cnap.reacBoxes(:,5) == i;
    cnap.reacBoxes(reacs_on_map,5) = cnap.nummaps+1;
    cnap.unsaved_changes = 1;
    cnap.local.boolGeneratedMapWasUpdated = 1;
    cnap = CNAaddMap(cnap,filename{i},['generated_map_' num2str(i-min(newmaps)+1,'%.2i')],dimensions(i),1);
end
disp('Placing reaction boxes');
for reacIndex = reacsToBePlacedOnNewMap(:)'
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


function str = getEQ( cnap, reacidx )
    zw=find(cnap.stoichMat(:,reacidx)<0);
    str = [strtrim(cnap.reacID(reacidx,:)) ': '];
    if(~isempty(zw))
            str=[str,  num2str(-cnap.stoichMat(zw(1),reacidx)),' ',deblank(cnap.specID(zw(1),:))];
            for j=2:length(zw)
                    str=[str,' + ',num2str(-cnap.stoichMat(zw(j),reacidx)),' ',deblank(cnap.specID(zw(j),:))];
            end
    else
        str = '';
    end
    if(~strcmp('mue',deblank(cnap.reacID(reacidx,:))))
        if cnap.reacMin(reacidx)>=0 && cnap.reacMax(reacidx)>0
            str=[str,' ==> '];
        elseif cnap.reacMin(reacidx)<0 && cnap.reacMax(reacidx)>0
            str=[str,' <=> '];
        elseif cnap.reacMin(reacidx)<0 && cnap.reacMax(reacidx)<=0
            str=[str,' <== '];
        end
    end
    zw=find(cnap.stoichMat(:,reacidx)>0);
    if(~isempty(zw))
            str=[str,num2str(cnap.stoichMat(zw(1),reacidx)),' ',deblank(cnap.specID(zw(1),:))];
            for j=2:length(zw)
                str=[str,' + ',num2str(cnap.stoichMat(zw(j),reacidx)),' ',deblank(cnap.specID(zw(j),:))];
            end                                                  
    end
end