function saveFigurePNG(folder_name_for_saving,figName)
%
% Created by Isabel Llorente-Garcia. March 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Save current figure as a .png file in a directory named "folder_name_for_saving",
% which already exists. The figure is saved with the name "figName".
% Both inputs are strings.
% Close the figure window at the end.
%
% Example: figName = strcat('frameAvg',image_label), with image_label='554'.

cd(folder_name_for_saving); % change to directory in which results are saved.
% Export the current figure window at screen size as a png into folder
% folder_name_for_saving.
set(gcf, 'PaperPositionMode', 'auto')  % Use screen size. (gcf=get current figure)
% h = get(gcf);
print('-dpng','-r300',figName)  % add -r300 (to save at 300 dpi, higher resolution) after -dpng to control resolution.
cd('..'); % go back to previous directory.
close; % deletes the current figure (many open figures take too much memory and make Matlab crash).
    