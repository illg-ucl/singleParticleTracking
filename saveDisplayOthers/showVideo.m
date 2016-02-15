function showVideo(image_label,start_frame,end_frame)
%
% Created by Isabel Llorente-Garcia, April 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Read image sequence (.sif, .dv, .tif or .mat data) and show as video.
%
% This function chooses an image sequence file which is in the current directory, 
% and shows a video (plot after plot) of frames from start_frame to end_frame.
% 
% example of how to call this function: showVideo('490',100,120) extracts and plots
% frames 100 to 120 of the chosen image sequence displaying them on a video.
%
%
% Inputs: 
% image_label: string that labels a given image sequence found in current
% folder (e.g.'513', '490', etc...) 
%


% % uigetfile opens a file dialog box to choose image data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% data_folder_path = strcat(path_data,file_data);
% disp(data_folder_path) % write .sif image path to command window.

% Read image sequence (.sif or .dv data)
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);

% To show a video:
% Loop through frames:
for p = start_frame:end_frame
    disp(['frame number: ',num2str(p)]) % print frame number to Command Window.
    imshow(image_data(p).frame_data,[],'Border','tight','InitialMagnification',150);
    pause(0.2);
end

