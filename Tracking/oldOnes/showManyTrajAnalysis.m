function processedManyTrajs = showManyTrajAnalysis(image_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo) 
%
% Isabel, Oct 2011.
% Analyse all trajectory data for a given image sequence (get mean square displacement (msd) vs Delta t plots and fits, and
% exponential decay of Intensity vs time, etc...), from trajectory
% 'n_traj_start' to 'n_traj_end'.
% Show results of trajectory analysis and save them, and also show the
% trajectory overlayed on the image sequence on a video to check.
% Save also results of trajectory analysis to an excel file (done within showTrajAnalysis.m).
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
% Eg.  T0 = showManyTrajAnalysis('101',5,7,5,0.04,35.333,1);  
% image '101' in folder, look at analysed trajectories 5 to 7, with frame 5
% being the time origin and 40ms between frames.
% Eg. to show only one trajectory (no. 8 of ATPase-GFP_101fullTrajs.xls, eg) do:
% T0 = showManyTrajAnalysis('101',8,8,5,0.04,35.333,1);
% To get vector with results labelled as "good trajectories" by the user
% do: [T0.IntensityAtTimeOrigin].*[T0.GoodTrajFlag]. (this does not work
% now).
% eg. analyse only traj number 2 for image 500 in 'cybD-mCherry-ATPase-GFp':  T500 = showManyTrajAnalysis('500',2,2,5,0.04,35.333,1);
% To get results, do T500{1}{i},   for i from 1 to 8.  
% ------------------------------------


%% PARAMETERS: 

% Set the "quickLook" parameter as 1 if you are just looking quickly at
% "good" trajectories and you don't want to be asked for user input as to
% whether the trajectory is "good" or not (to flag them), and you don't
% want to save the result structure as a .mat in a user specified folder
% either...
quickLook = 0;

%  PARAMETER minPointsTraj: Minimum number of data points that a trajectory must have in order to be
% analised:
% Note that this number needs to be at least 5 for all methods in "showTrajAnalysis.m" to work well.
minPointsTraj = 5;


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


%% Analyse all trajectory data (get msd): 
% (.xls file previously generated with functions 'FindTrajects' and 'linkTrajSegments'):
analysedAllTraj = analyseTraj(trajXlsPath,tsamp,minPointsTraj);
% 'analysedAllTraj' is a structure array with as many elements as
% analysed trajectories (the ones with at least 15 points in them), and
% with fields 'XLS', 'xvalues', 'yvalues', 'intensity', 'msd_unavg',
% 'timeabs', 'timerel', 'numel','minNumPointsInTraj', 'timediff', 'msd', 'errorMsd',
% 'errorMsdRelPercent' and 'disp'.

%'The total number of trajectories analysed (long enough) in the file is: 
n_trajs_analysed = length(analysedAllTraj);


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

n_good_tracking = 1; % initialise index for trajs with good tracking which are saved within loop.

for n = n_traj_start:n_traj_end
    
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
        % 'position' vector is [left, bottom, width, height].
        % left, bottom control the position at which the window appears when it pops.
        
        for k = 1:length(frames_list)
            
            frame = image_data(frames_list(k)).frame_data; % extract frame data which is stored in field 'frame_data'.
            frame = double(frame);
            
            imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
            hold on;
            
            plot(x_values(k),y_values(k),'o','Color','g','MarkerSize',12) % plot accepted spot centres in green.
            pause(0.1); % this pause is needed to give time for the plot to appear
            hold off;
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
        processedTraj = showTrajAnalysis(trajXlsPath,image_data,analysedAllTraj,n,start_frame,tsamp,pixelsize_nm);
        
        % output: structure array, starting at element 1, only save to results trajs with good tracking.
        processedManyTrajs{n_good_tracking} = processedTraj; % Function output: cell array, starting at element 1.
        
        % Note that saving is done within function "showTrajAnalysis.m". 
        
        n_good_tracking = n_good_tracking + 1;
    end
    
    
end


%% Save result (as .mat) in a folder specified by user.
% Do not save if we are just having a quick look at good trajectories
% (quickLook = 1):
if quickLook ~=1
    uisave('processedManyTrajs','procManyTraj')
end


