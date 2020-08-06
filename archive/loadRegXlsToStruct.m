function regulation = loadRegXlsToStruct(cnap,xlsFile)
    TableReads = loadSpecReacXLStoStrArray(xlsFile);
    
    [ridrow,ridcol,regSheet]   = ind2sub(size(TableReads),find(strcmp(strtrim(TableReads),'reg_ind')));
    [rdurow,rducol,~]   = find(strcmp(strtrim(TableReads(:,:,regSheet)),'reg_down_up'));
    [rbdrow,rbdcol,~]   = find(strcmp(strtrim(TableReads(:,:,regSheet)),'reg_bounds'));
    
    regulation.reg_ind = find(strcmp(cellstr(cnap.reacID),TableReads(ridrow+1,ridcol,regSheet)));
    
    regulation.reg_down_up(1,1) = str2double(cell2mat(TableReads(rdurow+1,rducol,regSheet)));
    regulation.reg_down_up(2,1) = str2double(cell2mat(TableReads(rdurow+2,rducol,regSheet)));
    
    lastrow = rbdrow+find(strcmp(strtrim(TableReads(rbdrow:end,rbdcol,regSheet)),''),1,'first')-1;
    % Else take last row with content
    if isempty(lastrow)
        lastrow = rbdrow+find(~strcmp(strtrim(TableReads(rbdrow:end,rbdcol,regSheet)),''),1,'last')-1;
    end
    for i = rbdrow+1:lastrow
        regulation.reg_bounds(i-1) = str2double(cell2mat(TableReads(i,rbdcol,regSheet)));
    end
end