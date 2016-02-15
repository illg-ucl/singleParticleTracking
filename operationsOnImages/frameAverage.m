function frame_avg = frameAverage(image_label,start_frame,end_frame,display_result,save_png)
% 
% Created by Isabel Llorente-Garcia, January 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function takes an image_label such as '513', '490', etc... 
% which corresponds to a certain .sif image sequence in current folder,
% and then calculates a frame average for frames between start_frame and
% end_frame.
% If end_frame input is 'end', it uses all frames.
% And shows result only if "display_result" is equal to 1.
% And saves resulting frame average as a png image only if input "save_png" is
% equal to 1.
%
% Example of how to call this function:
% A = frameAverage('470',5,500,0,0); uses
% image sequence 470 within the current folder
% and finds the frame average for frames 5 to 500.
% Another example: A = frameAverage('490',10,200,1,1); does the same for frames 10
% to 200, displays result too, and saves a .png of resulting avg image.
% -------------------------------------------------------


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------

if strcmp(end_frame,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   end_frame = numFrames;
end


%% Calculate frame average:

% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul = zeros(frame_Ysize,frame_Xsize);

for k = start_frame:end_frame  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul = frame_accumul + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg = frame_accumul/(end_frame-start_frame+1);

if display_result==1,
    imshow(frame_avg,[]); % display result for graphical inspection
end


%% Save/export frame-average image:

if save_png==1
    
    % Make new folder (new directory) to save frame-average result:
    new_folder_name = 'frameAverages';
    warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
    mkdir(new_folder_name); % make new directory.
    cd(new_folder_name); % change to that directory.
    
    % Export the current figure window at screen size as a png:
    % Result figure files are saved in the new folder within the same directory where the input .sif image was.
    figName = strcat('frameAvg',image_label); % name of figure file to save to.
    set(gcf, 'PaperPositionMode', 'auto')  % Use screen size.
    print('-dpng','-r300',figName)  % add -r300 (to save at 300 dpi, higher resolution) after -dpng to control resolution.
    
    cd('..'); % go back to previous directory.
    
    close; % deletes the current figure.
    
end
