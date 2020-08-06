function [table_reads sheets] = loadSpecReacXLStoStrArray( met_reac_set_path )
    % Returns a 3D string table that contains all information imported from an
    % xls-file. Dimentions are derived from worksheet with the biggest range
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 13, 2017
    filepath = which(which(met_reac_set_path));
    if verLessThan('matlab', '9.7')
        [~,sheets,~] = xlsfinfo(filepath);
        % Import data to 3D array [row,column,sheet]
        for i=1:length(sheets)
            [~,~,workSheet] = xlsread(filepath,char(sheets(i)));
            if ~isempty(workSheet)
                table_reads(1:size(workSheet,1),1:size(workSheet,2),i) = workSheet;
            end
        end
        table_reads = strTable(table_reads);
    else
        sheets = sheetnames(filepath);
        % Import data to 3D array [row,column,sheet]
        for i=1:length(sheets)
            workSheet = readcell(which(met_reac_set_path),'sheet',i);
            if ~isempty(workSheet)
                table_reads(1:size(workSheet,1),1:size(workSheet,2),i) = workSheet;
            end
        end
        table_reads = strTable(table_reads);
    end
end

function rawmat = strTable(rawmat)
    [rows, cols, shts] = size(rawmat);
    for row = 1:rows
        for col = 1:cols
            for sh = 1:shts
                if isnumeric(rawmat{row,col,sh})
                    if isnan(rawmat{row,col,sh}) | isempty(rawmat{row,col,sh})
                        rawmat(row,col,sh) = {''};
                    else
                        rawmat{row,col,sh} = num2str(rawmat{row,col,sh});
                    end
                elseif isa(rawmat{row,col,sh},'missing')
                    rawmat(row,col,sh) = {''};
                end
            end
        end
    end
end