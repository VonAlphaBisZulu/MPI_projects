function [compchoice, extra] = diag_example(cnap,options,diag_title, diag_text,checkbox_text)
% example usage:
% diag_example({'Choice 1'; 'Choice 2', 'Choice 3'},'Hello World','Say hello to the world','Do you want the extra?')

if nargin == 1
    diag_title = 'example';
    diag_text = 'Please make your choice';
    options   = {'Option 1'; 'Option 2'; 'Option 3'};
    checkbox_text = 'activate checkbox';
    cb = 1;
else
    if ~exist('diag_title','var')
        diag_title = 'Title';
    end
    if ~exist('diag_text','var')
        diag_title = 'Please make your choice';
    end
    if nargin < 4
        cb = 0;
    else
        cb = 1;
    end
end

compchoice = [];
extra = [];

%% Make Window
scrsize=get_screen_size();
windowsize = [400  140];
diag=figure('Units','pixels','NumberTitle','off','Name',diag_title,'Position',[scrsize(3)/2-150 scrsize(4)/2-120 windowsize ],'Menubar','none','NumberTitle','off');
text('String',diag_text,'Units','pixels','Position',[75, 120],'FontUnits','pixels','FontSize',18);
tb = uicontrol('Style','edit','Units','pixels','Parent',diag,'String', '###','HorizontalAlignment','left','Position',[20,80,360,20]);
% ========== Hier wird die autocomplete-Funktionalität hinzugefügt ============
set(tb,'UserData',struct('isAutoCompleted',0),'BusyAction','cancel','KeyPressFcn',{@autocompleteText, findjobj(tb,'persist') , evalin('base',[cnap.net_var_name '.reacID'])}) % if cnap is passed directly: cnap.reacID
% =============================================================================
set(diag,'Resize','off');
set(gca,'Visible','off','Position',[0,0,1,1]);
% Put in buttons and checkboxes
if cb
    cbox=uicontrol('Style','checkbox','Units','pixels','Position',[20,50,20,20],'Value',0);
    text('String',checkbox_text,'Units','pixels','Position',[45, 60],'FontUnits','pixels','FontSize',12,'ButtonDownFcn',@(varargin) set(cbox,'Value',~get(cbox,'Value')));
end
buttons = [];
for opt = 1:length(options)
    buttons = [buttons, uicontrol('Style','pushbutton','String',options(opt),'Units','pixels','FontUnits','pixels','FontSize',12,'UserData',[options(opt);num2cell(opt)],'callback',@execute_callback)];
end

% arrange buttons according to the available solvers
for button = 1:length(buttons)
    set(buttons(button),'Position',[10 + 380/length(buttons)*(button-1)+2, 10, (windowsize(1)-20-2*length(buttons))/length(buttons), 30])
end

waitfor(diag);

    function execute_callback(varargin)
        UData = get(varargin{1},'UserData');
        compchoice = UData{2};      % Get Choice
        if cb
            extra = logical(get(cbox,'Value'));        % Get checkbox
        end
        closereq;
    end
end