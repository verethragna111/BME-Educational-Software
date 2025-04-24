
currentDirectory = pwd;

subDirList = config_demo();

for it=1:length(subDirList)
    addpath(fullfile(currentDirectory,subDirList{it}));
end

demo_run()

%global data used;

%data = [];
%used = -1;

function demo_run()
    cols = [30,110,190,270,350];

    a = figure('Color',[0.8 0.8 0.8], ...
        'Colormap',jet, ...
        'Position',[155 150 320 130], ...
        'Tag','Fig1');

    % Popup menu
    popup = uicontrol('Parent', a, ...
        'Units','points',...
        'Style','popupmenu',...
        'Position', [cols(1) 70 180 20],...
        'String',  {'voltage clamp','current clamp'},...
        'Tag', 'PopUp1');

    % Pushbutton
    b = uicontrol('Parent',a, ...
        'Units','points', ...
        'Callback', {@run_gui, popup}, ... % Pass the popup handle to the callback
        'Position',[cols(1) 10 180 64], ...
        'FontSize', 20,...
        'String','run program', ...
        'Tag','Pushbutton1');
end

% Function to call GUI depending on the Popup menu selected
function run_gui(~, ~, popup)
    gui_selected = popup.Value; % Get the value of the popup menu

    if gui_selected == 1 
        voltageGUI();
    else 
        currentGUI();
    end 
end