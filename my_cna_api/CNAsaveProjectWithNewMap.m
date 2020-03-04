function cnap = CNAsaveProjectWithNewMap( cnap )
    %% saves map configuration changes to file (modifies app_para.m)
    % 
    % Usage:  CNAsaveProjectWithNewMap( cnap )
    %
    % Input: 
    %
    %   cnap:         CNA project object
    %
    % Output:
    %
    %   cnap:         CNA project object
    %
    % saves map configuration changes to file (modifies app_para.m). When
    % there were maps generated in the meantime, they are handled
    % accordingly. The old image file is replaced by the newly generated
    % one and the project fields are updated accordingly.
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 16, 2017
    
    res= questdlg(['Do you want to edit the projects app-para.m to store the added network map(s)? Path: '...
        ,cnap.path],'YES','YES','NO','NO');
    if strcmp(res, 'YES')
        %% =================== manage generated map(s) ========================
        % temporarily added generated maps are named "_generated.pcx". When
        % project is saved, the old "generated.pcx" file is replaced by the
        % newly generated one.
        if isfield(cnap,'local')
            if isfield(cnap.local,'boolGeneratedMapWasUpdated')
                if cnap.local.boolGeneratedMapWasUpdated
                    % Get path to current graph image file
                    [filepath,nameTempGenerated,ext] = fileparts(char(cnap.maps(strcmp(cnap.maps(:,1),'generated'),2)));
                    if exist([cnap.path,'/',filepath,'_generated',ext], 'file')==2
                        % Delete old generated graph
                        if exist([cnap.path,'/',filepath,'generated',ext], 'file')==2
                            delete(strcat(cnap.path,'/',filepath,'generated',ext));
                        end
                        % Copy newly generated graph to "replace" the old one
                        copyfile([cnap.path,'/',filepath,nameTempGenerated,ext],[cnap.path,'/',filepath,'generated',ext]);

                        % Update path to graph in project's map list
                        cnap.maps(strcmp(cnap.maps(:,1),'generated'),2) = {strcat(filepath,'generated',ext)};

                        % Delete temporary file
                        delete(strcat(cnap.path,'/',filepath,nameTempGenerated,ext));
                        cnap.local.boolGeneratedMapWasUpdated = 0;
                    end
                end
            end
        end
        % delete all generated maps if they're no longer used in project
        if ~any(strcmp(cnap.maps(:,1),'generated'))
            if exist([cnap.path,'/generated.pcx'],'file')
                delete([cnap.path,'/generated.pcx']);
            end
            if exist([cnap.path,'/generated.png'],'file')
                delete([cnap.path,'/generated.png']);
            end
            if exist([cnap.path,'/_generated.pcx'],'file')
                delete([cnap.path,'/_generated.pcx']);
            end
            if exist([cnap.path,'/_generated.png'],'file')
                delete([cnap.path,'/_generated.png']);
            end
        end

        % set flags
        cnap.local.boolGeneratedMapWasUpdated = 0;
        cnap.local.boolMapChangesInProject = 0;

        %% ==================== generate app_para.m ===========================
        % write everything to file
        % copied from neteditmaskeval, line 110ff (Jun 16, 2017) and adapted
        lfi=fopen([cnap.path '/app_para.m'],'w');
        if(lfi==-1)
            msgbox('Could not open file: app_para.m. Check permissions!');
            cd(cnan.cnapath);
            return;
        end
        fprintf(lfi,['epsilon=',num2str(cnap.epsilon),';\n']);
        fprintf(lfi,['basic_color=[',num2str(cnap.color1),'];\n']);
        fprintf(lfi,['cr_color=[',num2str(cnap.color2),'];\n']);
        fprintf(lfi,['br_color=[',num2str(cnap.color3),'];\n']);
        fprintf(lfi,['nbr_color=[',num2str(cnap.color1),'];\n']);
        fprintf(lfi,['text_color=[',num2str(cnap.textColor),'];\n']);
        fprintf(lfi,['macro_synth_color=[',num2str(cnap.macroSynthColor),'];\n']);
        fprintf(lfi,['macro_color=[',num2str(cnap.specBoxColor),'];\n']);
        fprintf(lfi,['box_reaction_width=[',num2str(cnap.reacBoxWidth),'];\n']);
        fprintf(lfi,['box_reaction_height=[',num2str(cnap.reacBoxHeight),'];\n']);
        fprintf(lfi,['box_macro_width=[',num2str(cnap.specBoxWidth),'];\n']);
        fprintf(lfi,['box_macro_height=[',num2str(cnap.specBoxHeight),'];\n']);
        fprintf(lfi,['fontsize_reaction=[',num2str(cnap.reacFontSize),'];\n']);
        fprintf(lfi,['fontsize_macro=[',num2str(cnap.specFontSize),'];\n']);
        fprintf(lfi,'fluxmaps={\n');
        for i=1:cnap.nummaps
            fprintf(lfi,['''',cnap.maps{i,1},''',''',strrep(cnap.maps{i,2},'\','\\'),'''\n']);
        end
        fprintf(lfi,'};\n');

        fclose(lfi);
        
    else
        disp('map configuration was not saved');
    end
end