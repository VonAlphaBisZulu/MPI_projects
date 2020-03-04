function deleteInvalidTextBoxes( cnap )
    % searches for text boxes, whose handles are not stored in the cna project
    % object, and deletes those text boxes. Affects only the gui.
    %
    % Philipp Schneider - schneiderp@mpi-magdeburg.mpg.de
    % Jun 16, 2017
    
    textBoxes = [];
    % get all matlab figures
    for fig = findall(cnap.figs(:,1),'Type','figure')'
        children = get(fig,'Children');
        for obj = children(:)'
            % get all UIControl objects / textboxes
            if isa(obj,'matlab.ui.control.UIControl')
                textBoxes = [textBoxes obj];
            end
        end
    end
    % if the textbox handles are not stored, delete the objects
    for tb = textBoxes
        if ~any(tb==cnap.reacBoxes(:,4))
            delete(tb);
        end
    end
end

