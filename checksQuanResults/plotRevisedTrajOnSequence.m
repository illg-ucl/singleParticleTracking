function plotRevisedTrajOnSequence(data_set_label,image_label,traj_numbers,pause_seconds)
% This function takes a data_set_label such as 'ATPase-GFP TIRF FRAP' or 'cybD-mCherry &ATPase-GFp', 
% an image_label such as '513', '490', etc... which corresponds to a certain .sif number of image sequence, 
% and certain trajectory numbers, traj_numbers, revises the centroids found
% for bright spots using function "findSpotCentre1frame" and then
% shows the image frames together with the revised trajectory results overlaid.
% The maximum trajectories that can be plotted is determined by the length of the list of colours, color_list, below.
%
% traj_numbers is a vector with the numbers of trajectories I want to plot. 
% Eg. traj_numbers = [1 2 5]. 
% Each trajectory contains lots of zeros as well as the relevant trajectory data.
%
% pause_seconds is the time in seconds between showing different frames 
%
% Several trajectories are plotted on the same figure: stuff..._image_label_paricle_trajectory_1, 
% stuff..._image_label_paricle_trajectory_2, etc.
%
% example of how to call this function:  
% plotRevisedTrajOnSequence('GFP-nuoF & mCherry-sdhC','470',[1:48],0.5) plots
% image sequence 470 within the folder 'GFP-nuoF & mCherry-sdhC', 
% and overlays trajectories 1 to 48 on top of the figure. Video is
% displayed with 0.5seconds in between each frame.
% Another example: plotRevisedTrajOnSequence('cybD-mCherry &ATPase-GFp','490',[16 17 18],0.5)

% define data (.sif image) and analysis (.xls trajectories) directories:
% dir_data = 'Z:\Leake\Heiko Data\';
dir_data = 'C:\Isabel\ExperimData\HeikoData\';
dir_analysis = 'C:\Isabel\DataAnalysis\oxphos\XueAnalysis\';
% print data and analysis directories on command window to guide user:
disp(' ') % empty line
disp(['The data directory (.sif images) is: ',dir_data])
disp(' ') % empty line
disp(['The analysis directory is: ',dir_analysis])

% uigetfile opens a file dialog box to choose data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% strcat('data (.sif image):','  ',path_data,file_data)
% open a file dialog box to choose analysis file:
% [file_analysis,path_analysis] = uigetfile({'*.xls'}, 'Chose analysis file (trajectory):');
% strcat('analysis file (.xls trajectory):','  ',path_analysis,file_analysis)

% list of colours for plotting each trajectory:
% colors are in this order: red, green, blue, yellow, pink, lighter blue,
% darker red, darker green, darker blue, mustard, darker pink, light blue,
% light orange, light greenish, darker green, salmon, another light blue, very light green,
% lighter grey, light gray, light yellow, light pink, light blue, orange
% the length of this list determines the total no. of trajectories we can plot (48 now).
color_list = [1 0 0;     0 1 0;   0 0 1;     1 1 0;     1 0 1;     0 1 1;...
               0.8 0 0; 0 0.8 0; 0 0 0.8; 0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8;...
               0.95 0.72 0.4; 0.86 0.95 0.6; 0.4 0.8 0.4; 0.9 0.6 0.5; 0.4 0.8 0.7; 0.85 0.92 0.99;...
               0.9 0.9 0.9; 0.7 0.7 0.7; 1 1 0.5; 1 0.5 1; 0.5 1 1; 1 0.6 0.2;...
               1 0 0;     0 1 0;   0 0 1;     1 1 0;     1 0 1;     0 1 1;...
               0.8 0 0; 0 0.8 0; 0 0 0.8; 0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8;...
               0.95 0.72 0.4; 0.86 0.95 0.6; 0.4 0.8 0.4; 0.9 0.6 0.5; 0.4 0.8 0.7; 0.85 0.92 0.99;...
               0.9 0.9 0.9; 0.7 0.7 0.7; 1 1 0.5; 1 0.5 1; 0.5 1 1; 1 0.6 0.2];

% error control:
if size(traj_numbers,2)>length(color_list)
    disp(' ') % empty line
    disp(['The max no. of trajectories you can plot is: ',num2str(length(color_list))])
    disp(' ') % empty line
    error('plotTraj:tooMany','you are trying to plot too many trajectories...')
end
           
% Change to the DataAnalysis directory:
cd(dir_analysis);

analysis_folder = strcat(cd,'\',data_set_label);
cd(analysis_folder); % change to the data folder
TIRF_analysis_folder = dir('*TIRF');  % returns a structure with all directories with the word TIRF on their name. The structure has fields 'name', 'date', 'bytes', etc. 
cd(TIRF_analysis_folder.name); % The directory name is stored in the 'name' field of the structure.

pre1 = strcat('*',image_label);
image_folder_struct = dir(pre1); % image_folder_struct is also a structure with field 'name'.
image_folder_path = strcat(cd,'\',image_folder_struct.name);

% Get the trajectories:
gral_traj_path = strcat(image_folder_path,'\','intensity\'); % general trajectories path
cd(gral_traj_path);

% select trajectory data:
all_traj_data = cell(1,size(traj_numbers,2));   % creates cell array of empty matrices to store trajectories. 
% No. of columns in array is the number of trajectories, i.e.,
% size(traj_numbers,2).
all_frame_numbers = [];
% now loop through trajectories:
for n = 1:size(traj_numbers,2)
    pre2 = strcat('*_',num2str(traj_numbers(n)),'.xls'); % pre2 is a string to look for in order to list files in the directory containing that string
    traj_name_struct = dir(pre2); % this is an structure with the field 'name' containing the name of the trajectory file
    traj_path = strcat(cd,'\',traj_name_struct.name); % string, name of excel trajectory file
    traj_data = xlsread(traj_path); % to open excel files need to use xlsread
    selected_rows = find(traj_data(:,1)~=0); % eliminate rows of zeros (with zero in first column).
    selected_traj_data = traj_data(selected_rows,[1,3,4,14,15]); % select meaningful columns: frame no. and x-y data:
    % In the traj_data, the column 1 is the frame number, columns 3 and 4
    % are global x and y on the top channel (red), and columns 14 and 15 are
    % global x and y on the bottom channel (green).
    all_traj_data(1,n) = {selected_traj_data}; % asign data from trajectory n to element n of cell array "all_traj_data".
    all_frame_numbers = [all_frame_numbers;all_traj_data{n}(:,1)]; % all frame numbers from trajectories, together in a vector.
end

start_frame = min(all_frame_numbers); % the function only plots from the first frame with a trajectory to the last one.
end_frame = max(all_frame_numbers);
disp(' ') % empty line
disp(['The start frame for the plot will be ',num2str(start_frame)])
disp(['The end frame for the plot will be ',num2str(end_frame)])

% % help to debug: to get the selected numeric data for each trajectory:
% all_traj_data{1}
% all_traj_data{2}
% cell2mat(all_traj_data(3))


% read in the .sif image data:
% Change to the Data directory:
cd(dir_data);

data_folder = strcat(cd,'\',data_set_label);
cd(data_folder); % change to the data folder
TIRF_data_folder = dir('*TIRF');  % returns a structure with all directories with the word TIRF on their name. The structure has fields 'name', 'date', 'bytes', etc. 
cd(TIRF_data_folder.name); % The directory name is stored in the 'name' field of the structure.

pre3 = strcat('*_',image_label,'.sif'); % pre3 is a string to look for in order to find files in the directory containing that string
data_folder_struct = dir(pre3); % data_folder_struct is also a structure with field 'name'.
data_folder_path = strcat(cd,'\',data_folder_struct.name);
disp(' ') % empty line
disp(data_folder_path) % write .sif image path to command window.

% first get size of .sif image file: (see IO_Input folder, SifFunctions.txt, or page 95 my notebook 1):
[ReturnCode, numFrames, ImageSize, TotalAcquisitionSize]=GetAndorSifSize(data_folder_path,0);
disp(' ') % empty line
disp(['The total number of frames in this image sequence is: ',num2str(numFrames)]) %output total number of frames to command window.
% numFrames is the length of the sequence, the number of frames.
% ImageSize is the size of the image, e.g. 512*512.
% TotalAcquisitionSize is numFrames*ImageSize

% read .sif image data. sifData is an array of 1x1 structures, as many
% columns as frames on image sequence. It reads frames 1 to numFrames:
[sifData] = read_sif_data_direct(data_folder_path,numFrames,ImageSize,1,numFrames);


% loop through frames plotting frame and corresponding trajectories for that frame. 
% Makes a video and only plots between the first frame with a trajectory and the last one:
for k = start_frame:end_frame
    % to plot all frames do instead: for k = 1:length(sifData)
    
p = sifData(k); % each column of sifData is a 1x1 structure with the data saved in the field 'sliceData'.
frame = p{1}.sliceData; % extract frame data which is saved in the field 'sliceData'.
framebis = mat2gray(frame); % transform matrix into image to do image operations on it.
framebisbis = imrotate(framebis,90); % strange..., it needs to be rotated by 90 degrees...
imshow(framebisbis,[],'Border','tight','InitialMagnification',150); % rotate image by 90 degrees and show it scaled between its min and max values ([]).
hold;
disp(['frame number: ',num2str(k)]) % print frame number to Command Window.

    % for each frame, loop throuh all trajectories:
    for m = 1:size(traj_numbers,2)
        
    % plot trajectory data in a scatter plot on top of the corresponding image frame:
    % all_traj_data is a cell array containing all trajectories.
    % all_traj_data{m} is the mth trajectory. Its first column contains the frame number
    select_traj_row = find(all_traj_data{m}(:,1)==k); % find row of trajectory data values that corresponds to frame k.
        if length(select_traj_row)==1
             top_x_data = all_traj_data{m}(select_traj_row,2); % second column is x_top data
             top_y_data = all_traj_data{m}(select_traj_row,3); % third column is y_top data
             bottom_x_data = all_traj_data{m}(select_traj_row,4); % fourth column is x_bottom data
             bottom_y_data = all_traj_data{m}(select_traj_row,5); % fourth column is y_bottom data
             % Now revise trajectories using function findSpotCentre1frame:
             % use Quan's trajectory results (spot centre estimates) as
             % initial estimates and then iterate to find spot centre.
             % Image subarray ROI is a square of size 17x17 pixels (halfwidth is
             % 8 pixels), inner circular mask that moves inside the fixed 17x17
             % square has a radius of 5 pixels and the applied gaussian
             % mask has a sigma of 3 pixels:
             [top_x_revis top_y_revis] = findSpotCentre1frame(framebisbis,top_x_data,top_y_data,8,5,3);
             [bottom_x_revis bottom_y_revis] = findSpotCentre1frame(framebisbis,bottom_x_data,bottom_y_data,8,5,3);
             % Plot revised results:
             plot(top_x_revis,top_y_revis,'o','Color',color_list(m,:),'MarkerSize',10);
             plot(bottom_x_revis,bottom_y_revis,'o','Color',color_list(m,:),'MarkerSize',10);
        end     
    end
    hold off;
    
pause(pause_seconds)

    
end