function saveTrackVideos(image_label,n_traj_start,n_traj_end,minPointsTraj,saveVideo_noCircle) 
%
% Created by Isabel Llorente-Garcia, August 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Plot and save video (.avi file) of image with tracks overlaid on top (or with no tracks overlaid) for a given image sequence 'image_label', from trajectory
% 'n_traj_start' to 'n_traj_end'.
% 
% It works well for short tracks. (It is similar to showManyTrajAnalysis2.m.)
%
% NOTE: before running this function you should move into a directory which
% contains both the .sif/video image (labelled by 'image_label') which has previously been analysed with
% 'FindTrajects.m' and 'linkTrajSegments.m' to produce the .xls file which
% contains the trajectory results, which should be in the same directory as the .sif file.
%
% INPUTS: 
% - 'image_label' string that labels the image sequence under analysis, e.g. '101'.
% - 'n_traj_start': first trajectory we want to analyse and check.
% - 'n_traj_end': last trajectory we want to analyse and check. If the
% string 'end' is entered, we go through to the last analysed trajectory.
% - minPointsTraj: minimum number of data points that a trajectory must have in order to be
% analised. A value of 3 is used for "showTrajAnalysis2.m" and
% "showManyTrajAnalysis2.m".  % Use last input minPointsTraj=6 for
% Oxphos data.
% A value of at least 6 needs to be used for all methods in
% "showTrajAnalysis.m" (and therefore "showManyTrajAnalysis.m") to work
% well.
% - saveVideo_noCircle: 1 if you also want to save a video of the frames
% with no circle overlaid on top of bright spot. Input 0 if you don't want
% to also save that video.
% -----------------
% IMPORTANT NOTE!!!: The values of all inputs: "image_label",
% "n_traj_start", "n_traj_end" and "minPointsTraj" here need to be the same
% as those used later on for functions "showManyTrajAnalysis.m" or "showManyTrajAnalysis2.m", 
% otherwise the trajectory numbers will be different!.
% -----------------
%
% OUTPUT: 
% Video output files .avi are saved to output folder, either with overlaid
% track or without it. One/two video/s per track.
%
% Example of how to call this function:
% saveTrackVideos('533',1,'end',6,0); 
% Use last input minPointsTraj=6 for Oxphos data.
% ------------------------------------

% PARAMETERS:
% Save also a video of frames with no circle overlaid on top of image:
saveVideo_noCircle = 0;
% saveVideo_noCircle = 1;

%% Get path for trajectory data (excel file):

% You need to be in the correct directory before running the function!!!!
% Find paths in current folder which contain 'image_label' string:
trajXlsPath0 = dir(strcat('*',image_label,'*.xls')); % Trajectory data path (excel file with the full trajectories as returned by function "linkTrajSegments.m").
% Error control:
if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
    error('Check you are in the correct directory and run again. No .xls file found for that image number. Make sure image number is in between quotes ''.'); 
end
trajXlsPath = trajXlsPath0.name;
% trajXlsPath0 is a structure and the file names is stored in the field 'name'.


%% Create output directory for saving video (.avi) results:

% Make new folder (new directory) to save trajectory analysis results:
initial_folder_path = cd;
% Create new folder for outputs inside current folder:
output_folder_name = strcat('videosAvi_allTrajs_',image_label);
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(output_folder_name); % make new directory.
output_folder_path = strcat(initial_folder_path,'\',output_folder_name);


%% Analyse all trajectory data (get msd): 
% (.xls file previously generated with functions 'FindTrajects' and 'linkTrajSegments'):
analysedAllTraj = analyseTraj(trajXlsPath,1,minPointsTraj);
% I have used tsamp = 1 here, so the time is in units of frames (irrelevant since this function is only used for a visual check). 
% 'analysedAllTraj' is a structure array with as many elements as
% analysed trajectories, and
% with fields 'XLS', 'xvalues', 'yvalues', 'intensity', 'msd_unavg',
% 'timeabs', 'timerel', 'numel','minNumPointsInTraj', 'timediff', 'msd', 'errorMsd',
% 'errorMsdRelPercent' and 'disp'.

if isempty(analysedAllTraj)
    disp('The total number of long enough trajectories (to analyse) for this file is zero.');
    disp('Exiting program');
    return % exits function.
end

%'The total number of trajectories analysed (long enough) in the file is: 
n_trajs_analysed = length(analysedAllTraj);
disp(['The number of long enough tracks to go through for this image sequence is: ' num2str(n_trajs_analysed)])

if strcmp(n_traj_end,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   n_traj_end = n_trajs_analysed;
end


%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.



%% Loop selected trajectories:

good_track_numbers = []; % initialise empty vector where I will store numbers of tracks labelled as "good" by the user.

clear mex % close all open .avi files (this avoids errors).

% Move into output folder to save videos:
cd(output_folder_path);

for n = n_traj_start:n_traj_end

    % n

    % Close any pre-existing figures:
    close(findobj('Tag','Trajectory results'));
    % close(all); % close all figures.

    % Show video of trajectory overlaid on actual image:
    frames_list = analysedAllTraj(n).frame; % list of frame numbers in trajectory n.
    x_values = analysedAllTraj(n).xvalues0; % list of original x centres of spots in trajectory n.
    y_values = analysedAllTraj(n).yvalues0; % list of original y centres of spots in trajectory n.


    % Loop through frames in each trajectory analysed:
    fig = figure('Tag','Data video','units','inches','position',[12 4 6 6]); % Figure number 2.
    % 'position' vector is [left, bottom, width, height].
    % left, bottom control the position at which the window appears when it
    % pops.
    set(fig,'DoubleBuffer','on');
    set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')
    
    video_filename = strcat('video',image_label,'_traj',num2str(n),'.avi'); % filename of video with overlaid circle on top of bright spot.
    % Create and open avi file for saving frames onto it later:
    mov = avifile(video_filename); % this saves the .avi file in the current directory automatically. Movie with circle overlaid onto detected spot.
    
    if saveVideo_noCircle == 1
        % Save also a video without overlaid circle:
        video_filename_noCircle = strcat('video',image_label,'_traj',num2str(n),'_noCircle','.avi');
        % Create and open avi file for saving frames onto it later:
        mov_noCircle = avifile(video_filename_noCircle); % this saves the .avi file in the current directory automatically. Movie with no circle overlaid.
    end
    
    
    for k = 1:length(frames_list)

        frame = image_data(frames_list(k)).frame_data; % extract frame data which is stored in field 'frame_data'.
        frame = double(frame);
        
        imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
        
        if saveVideo_noCircle == 1
            % Save frame with no circle overlaid on top of bright spot, to video mov_noCircle:
            F_noCircle = getframe(gca);
            mov_noCircle = addframe(mov_noCircle,F_noCircle);
        end
        
        hold on;
        plot(x_values(k),y_values(k),'o','Color','g','MarkerSize',12) % plot accepted spot centres in green.
        pause(0.3); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
        hold off;

        % Save frame with circle overlaid, to the other video, mov:
        F = getframe(gca);
        mov = addframe(mov,F);
        
    end
    
    % Close avi files:
    mov = close(mov); 
    if saveVideo_noCircle == 1
        mov_noCircle = close(mov_noCircle);
    end
    
    disp(['Track number ',num2str(n),' out of ' num2str(n_trajs_analysed) ' tracks:'])
    
    % Flag saying if trajectory is a good one or not (bgnd point, not good
    % tracking, etc.).
    % Request user input if you want a pause in between videos of different tracks:
%     good_tracking_flag = input('Press any key to continue...'); 
       
    close(findobj('Tag','Data video')); % close video figure;
    
end

% Return to previous folder:
cd(initial_folder_path);

% --------------------