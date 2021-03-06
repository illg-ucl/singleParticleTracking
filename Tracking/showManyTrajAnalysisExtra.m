function processedManyTrajs = showManyTrajAnalysisExtra(image_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo,minPointsTraj,n_ExtraFrames) 
%
% Created by Isabel Llorente-Garcia, Oct 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function is similar to showManyTrajAnalysis2.m, but the "Extra"
% means that a few (n_ExtraFrames) extra points are added at the end of each track to try
% and capture the background level after a bright spot has dissappeared.
%
% Analyse all trajectory data for a given image sequence (intensity vs time, mean square displacement (msd) vs Delta t plots and fits, and
% exponential decay of Intensity vs time, etc...), from trajectory
% 'n_traj_start' to 'n_traj_end', but only for the trajectory numbers selected as
% "good" ones by eye, i.e., those within files
% "good_track_nums_image_label.mat". These files are generated by function
% goThroughTracksVideo(image_label,n_traj_start,n_traj_end,minPointsTraj).
% So it makes use of a list of "good track" numbers contained in a .mat file in the current directory (see below).
%
% Show results of trajectory analysis and save them, and also show the
% trajectory overlayed on the image sequence on a video to check.
% Save also results of trajectory analysis to an excel file (done within showTrajAnalysisExtra.m).
% 
% This works for short tracks (as opposed to showManyTrajAnalysis.m and
% showTrajAnalysis.m, which only work for tracks with at least 6 points in
% them). This function uses showTrajAnalysisExtra.m.
%
% NOTE: before running this function you should move into a directory which
% contains both the .sif image (labelled by 'image_label') which has previously been analysed with
% 'FindTrajects.m' and 'linkTrajSegments.m' to produce the .xls file which
% contains the trajectory results, which should be in the same directory as the .sif file.
%
% INPUTS: 
% - 'image_label' string that labels the image sequence under analysis, e.g. '101'.
% - 'n_traj_start': first trajectory we want to analyse and check.
% - 'n_traj_end': last trajectory we want to analyse and check. If the string 'end' is entered, we go through to the last analysed trajectory.
% - 'start_frame' is the number of frame considered as the origin of time (as t=0), in frames. It is the first frame for which
% the shutter is fully open and we detect fluorescence (written in my notebook when I analysed each image sequence).
% - 'tsamp' is the sampling time, used to calibrate the absolute time, to go from frames to time in seconds. 
% It is the time between frames in seconds. Use tsamp = 1 for the time to be in units of frames. A proper calibration
% have tsamp = 40*10^(-3), i.e., 40ms per frame, for example.
% start_frame*tsamp is therefore the absolute time origin in seconds.
% - 'pixelsize_nm': pixel size in nm (35.333nm for OXPHOS data).
% - 'showVideo' is an input parameter to show a video of the trajectory
% overlaid on the image sequence or not.
% - minPointsTraj: Minimum number of data points that a trajectory must have in order to be
% analised (default minPointsTraj = 3).
% Note that this number needs to be at least 3 for all methods in "showTrajAnalysisExtra.m" to work well.
% - n_ExtraFrames: number of extra frames added at the end of each track to
% try and capture the background level.
%
% OUTPUT: 'processedManyTrajs' is a cell array with as many elements (processedManyTrajs{n}) as
% analysed trajectories (those with good tracking only).
% The first element within trajectory {n} in the cell array, {n}{1}, is a structure with fields:
% fieldnames(processedManyTrajs{1}{1}):
%     'good_tracking_flag'
%     'XLSpath'
%     'minNumPointsInTraj'
%     'ImageSizeHorizPix'
%     'ImageSizeVerticPix'
%     'AnalysedTrajNum'
%     'OriginalTrajNum'
%     'FirstTrajFrame'
%     'LastTrajFrame'
%     'NumDataPoints'
%     'TrajStartTime'
%     'TrajEndTime'
%     'TrajDuration'
%     'TimeBetweenFrames'
%     'FrameForTimeOrigin'
%     'AbsTimeOrigin'
%     'pixelsize_nm'
%     'Track_meanX_0'
%     'Track_meanY_0'
%     'Track_meanX_cellOrigin'
%     'Track_meanY_cellOrigin'
%     'Track_meanXvalue_cellCoords'
%     'Track_meanYvalue_cellCoords'
%     'TopOrBottom'
%     'WindowWidthCKfilter'
%     'WeightExponentCKfilter'
%     'useFiltered'
%     'numBins'
%
% fieldnames(processedManyTrajs{1}{2}):
%     'frame'
%     'timeabs'
%     'intensity'
%     'xvalues0'
%     'yvalues0'
%     'xvalues'
%     'yvalues'
%     'x_values_cell'
%     'y_values_cell'
%
% fieldnames(processedManyTrajs{1}{3}):
%     'I0_fit'
%     'stDev_I0'
%     'tau_fit'
%     'stDev_tau'
%     'rsq_fit_I'
%     'I0_fit_wo'
%     'stDev_I0_wo'
%     'Ioffset_fit_wo'
%     'stDev_Ioffset_wo'
%     'tau_fit_wo'
%     'stDev_tau_wo'
%     'rsq_fit_I_wo'
%     'I_line_offset'
%     'stDev_I_line_offset'
%     'I_line__slope'
%     'stDev_I_line__slope'
%     'I_line_rsq'
%     'Npoints_for_I_line_fit'
%
% fieldnames(processedManyTrajs{1}{4}):
%     'timediff'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% fieldnames(processedManyTrajs{1}{5}):
%     'BinCentresPos'
%     'FreqCountsPos'
%
% fieldnames(processedManyTrajs{1}{6}):
%     'IstepAxis'
%     'PowerSpectrum'
%
% fieldnames(processedManyTrajs{1}{7}):
%     'Isteps'
%     'power_Fourier'
%
% fieldnames(processedManyTrajs{1}{8}):
%     'xsize_default_cell_region'
%     'ysize_default_cell_region'
%     'local_cell_region_xleft'
%     'local_cell_region_xright'
%     'local_cell_region_ytop'
%     'local_cell_region_ybottom'
%     'centre_of_mass_X'
%     'centre_of_mass_Y'
%     'coord_transform_angle'
%     'cell_width'
%     'cell_length'
%     'traj_in_pole'
%
% Eg.  T0 = showManyTrajAnalysis2('101',5,7,5,0.04,35.333,1);  
% image '101' in folder, look at analysed trajectories 5 to 7, with frame 5
% being the time origin and 40ms between frames.
% Eg. to show only one trajectory (no. 8 of ATPase-GFP_101fullTrajs.xls, eg) do:
% T0 = showManyTrajAnalysis2('101',8,8,5,0.04,35.333,1);
% To get vector with results labelled as "good trajectories" by the user
% do: [T0.IntensityAtTimeOrigin].*[T0.GoodTrajFlag]. (this does not work
% now).
% eg. analyse only traj number 2 for image 500 in 'cybD-mCherry-ATPase-GFp':  T500 = showManyTrajAnalysis2('500',2,2,5,0.04,35.333,1);
% To get results, do T500{1}{i},   for i from 1 to 8.  
% ------------------------------------


%% PARAMETERS: 

% Set the "quickLook" parameter as 1 if you are just looking quickly at
% "good" trajectories and you don't want to be asked for user input as to
% whether the trajectory is "good" or not (to flag them), and you don't
% want to save the result structure as a .mat in a user specified folder
% either...
quickLook = 0;


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


%% Create new directory for saving trajectory-analysis result structure for this image sequence:

% Make new folder (new directory) to save trajectory analysis results:
pos1 = strfind(trajXlsPath,'fullTrajs.xls'); % position of the start of the string 'fullTraj.xls' in the xls input file name.
new_folder_name = strcat(trajXlsPath(1:(pos1-1)),'_Extra'); % Take the name of the input excel file (with the end bit 'fullTraj.xls' removed) as the new folder name.
% Note that the directory "new_folder_name" is created by function
% showTrajAnalysisExtra.m when called from this function.


%% Get path for .mat file with the numbers of the good tracks:

% A file of "good" track numbers with a name "good_track_nums_1757.mat" (eg) is produced by
% function "goThroughTracksVideo.m":
good_tracks_path = dir(strcat('*good_track_nums','*',image_label,'*.mat'));
load(good_tracks_path.name); % This loads the structure "good_tracks" onto the workspace.
% good_tracks is a structure with fields:
% good_tracks.image_label = image_label; input value.
% good_tracks.n_traj_start = n_traj_start; input value.
% good_tracks.n_traj_end = n_traj_end; input value.
% good_tracks.minPointsTraj = minPointsTraj, input value.
% good_tracks.track_numbers: is a row vector containing the numbers of tracks considered as
% "good" by the user after seing the videos of the track overlaid on the
% image sequence.


%% Analyse all trajectory data (get msd): 
% (.xls file previously generated with functions 'FindTrajects' and 'linkTrajSegments'):
analysedAllTraj = analyseTraj(trajXlsPath,tsamp,minPointsTraj);
% 'analysedAllTraj' is a structure array with as many elements as
% analysed trajectories, and
% with fields 'XLS', 'xvalues', 'yvalues', 'intensity', 'msd_unavg',
% 'timeabs', 'timerel', 'numel','minNumPointsInTraj', 'timediff', 'msd', 'errorMsd',
% 'errorMsdRelPercent', 'disp', 'SNR','BgNoiseStd','IbgAvg','IinnerTot','rsqFit'.

if isempty(analysedAllTraj)
    disp('The total number of long enough trajectories (to analyse) for this file is zero.');
    disp('Exiting program');
    return % exits function.
end

%'The total number of trajectories analysed (long enough) in the file is: 
n_trajs_analysed = length(analysedAllTraj);


if strcmp(n_traj_end,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   n_traj_end = n_trajs_analysed;
end


%% Error control: IMPORTANT!!
% Check that the "minPointsTraj" (min no. of points in track for it to be
% analysed) and all other inputs for function "goThroughTracksVideo.m" were
% the same as for this function:
if (good_tracks.minPointsTraj ~= minPointsTraj) || ...
        (strcmp(good_tracks.image_label,image_label)~=1) || ...
        (good_tracks.n_traj_start ~= n_traj_start) || ...
        (good_tracks.n_traj_end ~= n_traj_end)
    disp('ERROR: parameters minPointsTraj, image_label, n_traj_start, n_traj_end must be the same as those used to create list of "good track" numbers with function goThroughTracksVideo.m or goThroughTracksVideo2.m. Exiting function...')
    return % exit function.
end
% row vector with  numbers of "good" tracks:
good_track_numbers = good_tracks.good_track_numbers;

disp(strcat('The number of "good" tracks in the good_track_nums file is:',num2str(length(good_track_numbers))))


%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.

% image_data is then an input to function showTrajAnalysisExtra.m.


%% Loop selected trajectories:

n_good_tracking = 1; % initialise index for trajs with good tracking which are saved within loop.

for n = n_traj_start:n_traj_end
    
    % Check if track number n is one of the "good" ones:
    B = ismember(good_track_numbers,n); % result is a vector with zeros at all positions except at the position of n in vector good_track_numbers, if it is a "good" one.
    % sum(B) is equal to 0 if "n" is not a "good" track, and equal to 1 if
    % "n" is on the list of "good" track numbers.
    
    if sum(B) == 1 % If track number "n" is a "good" one:
        
        % Close any pre-existing figures:
        close(findobj('Tag','Trajectory results'));
        % close(all); % close all figures.
        
        % Show video of trajectory overlaid on actual image:
        frames_list = analysedAllTraj(n).frame; % list of frame numbers in trajectory n.
        x_values = analysedAllTraj(n).xvalues0; % list of original x centres of spots in trajectory n.
        y_values = analysedAllTraj(n).yvalues0; % list of original y centres of spots in trajectory n.
        
        % Show video of trajectory overlaid on actual image:
        if showVideo == 1 % 'showVideo' is input parameter.
            
            % Loop through frames in each trajectory analysed:
            figure('Tag','Data video','units','inches','position',[12 4 6 6]); % Figure number 2.
            
            %         % Create video file:
            %         set(gca,'nextplot','replacechildren');
            
            % 'position' vector is [left, bottom, width, height].
            % left, bottom control the position at which the window appears when it pops.
            
            % Show the original track:
            for k = 1:length(frames_list)
                
                frame = image_data(frames_list(k)).frame_data; % extract frame data which is stored in field 'frame_data'.
                frame = double(frame);
                
                imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
                hold on;
                
                plot(x_values(k),y_values(k),'o','Color','g','MarkerSize',12) % plot accepted spot centres in green.
                pause(0.3); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
                hold off;   
            end
            
            % Show the added extra frames
            for kk = 1:n_ExtraFrames
                
                if (frames_list(k)+kk) <= numFrames % Error control, make sure we don't try to extend the track past the last frame.
                    
                    frame = image_data(frames_list(k)+kk).frame_data; % extract frame data which is stored in field 'frame_data'.
                    frame = double(frame);
                    
                    imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
                    hold on;
                    
                    plot(x_values(end),y_values(end),'o','Color','g','MarkerSize',12) % plot green circle at last tracked position.
                    pause(0.3); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
                    hold off;
                    
                end
            end
            
        end
        
        % For a quick analysis of "good" trajectories, quickLook = 1, no user
        % input requested and trajectories flagged as "good tracking":
        if quickLook ==1
            good_tracking_flag = 1;
        else
            % CHECK: skip the following visual check of every track or not.
            good_tracking_flag = 1;
            % good_tracking_flag = input('Is the tracking "good" for this trajectory? (1 for "yes", anything else for "no"): '); % request user input.
            % flag saying if trajectory is a good one or not (bgnd point, not good tracking, etc).
        end
        
        close(findobj('Tag','Data video')); % close video figure;
        
        if good_tracking_flag == 1
            % only analyse n-th trajectory if tracking is good (and folder created for good-tracking trajectories only).
            % "good_tracking_flag" is added to result structure in
            % "showTrajAnalysis.m", always 1, because if tracking is no good, the
            % trajectory is not analysed.
            
            % Analyse n-th trajectory data, produce result plots and save them:
            processedTraj = showTrajAnalysisExtra(trajXlsPath,image_data,analysedAllTraj,n,start_frame,tsamp,pixelsize_nm,n_ExtraFrames);
            
            % output: structure array, starting at element 1, only save to results trajs with good tracking.
            processedManyTrajs{n_good_tracking} = processedTraj; % Function output: cell array, starting at element 1.
            
            % Note that saving is done within function "showTrajAnalysis.m".
            
            n_good_tracking = n_good_tracking + 1;
        end
        
    end
end



%% Save result (as .mat) in output folder which contains track results for this image sequence:

cd(new_folder_name) % move into folder corresponding to track results for this image sequence.
output_filename = strcat('procManyTraj_Extra_',image_label);
save(output_filename,'processedManyTrajs'); % save variable processedManyTrajs to a .mat file
cd('..') % go back to previous folder.

% Save result (as .mat) in a folder specified by user input:
% % Do not save if we are just having a quick look at good trajectories
% % (quickLook = 1):
% if quickLook ~=1
%     uisave('processedManyTrajs','procManyTraj')
% end