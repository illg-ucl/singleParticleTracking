function tracks = analyseTraj(file,tsamp,minPointsTraj)
% Analyse all trajectory data and calculate mean square displacements.
%
% Function created by Alex Oct 2011 (formerly called 'parseXLS2mat.m').
% Modified and commented by Isabel October-2011.
% If you use this code please acknowledge Alex Robson and Isabel Llorente-Garcia in your
% publications.
%
%
% Inputs:
% 'file' is the path to an excel file with the full trajectories as returned by
% function "linkTrajSegments.m".
% 'tsamp' is the sampling time, used to calibrate the absolute time, to go from frames to time in seconds. 
% It is the time between frames in seconds. 
% 'minPointsTraj' : minimum number of data points that a trajectory must have in order to be
% analised.
%
% The output 'tracks' is a structure array with as many elements as
% analysed trajectories (the ones with at least minPointsTraj in them), and
% with fields 'trajNumber','xvalues0', 'yvalues0', 'xvalues', 'yvalues', 'intensity', 'I0_IspotFit', 'sigma_IspotFit', 'msd_unavg',
% 'frame', 'timeabs', 'timerel', 'numel', 'deltaTime', 'msd', 'errorMsd',
% 'errorMsdRelPercent' and 'disp' (see end of file for more details).
%
% Example of how to call this function:
% pathTraj = C:\Documents and Settings\llorentegarcia\Desktop\egTrajs.xls;
% using time in frames (tsamp = 1). (A proper calibration have tsamp = 40*10^(-3), i.e., 40ms per frame, for example.)
% t1 = analyseTraj(pathTraj,1); 
% or:
% pathTraj = 'C:\Isabel\ExperimData\HeikoData\ATPase-GFP\ATPase-GFP TIRF\ATPase-GFP_89fullTrajs.xls';
% sampling_time = 0.04; 
% t1 = analyseTraj(pathTraj,sampling_time);


% Dependencies: xlsread.mat
% Check dependencies:
if 2~=exist('xlsread') 
    error('Check file dependencies - you need to install xlsread'); 
end

% Open and read the data:
[NUMERIC,TXT,RAW]=xlsread(file,'Track results'); % import the data in the sheet named 'Track results'.

% Import the column heads and assign ID
colheads = TXT;

% The column titles are 'CentreX','CentreY','IspTot','rsqFit','SigmaFit','I0Fit','BgNoiseStd','IbgAvg','IbgTot'
% 'SNR','IinnerTot','ClipFlag','noConverge','TrajNumber','FrameNumber','SpotNumber'.

% Generate ID: ID is a structure with fiels with the same names as the
% column titles, and each has an ID value of 1, 2, 3, etc (see below).
for i=1:numel(colheads) 
    ID.(colheads{i}) = find(strcmp(TXT,colheads{i})); 
end
% eg. ID = 
%                             CentreX: 1
%                             CentreY: 2
%                              IspTot: 3
%                              rsqFit: 4
%                            SigmaFit: 5
%                               I0Fit: 6
%     bg_noise_offset_afterBGsubtract: 7
%                          BgNoiseStd: 8
%                              IbgAvg: 9
%                              IbgTot: 10
%                                 SNR: 11
%                           IinnerTot: 12
%                            ClipFlag: 13
%                          noConverge: 14
%                          TrajNumber: 15
%                         FrameNumber: 16
%                          SpotNumber: 17

% Write into a structure file as used by Alex

% The trajectory column:
traj = NUMERIC(:,ID.TrajNumber); % NUMERIC is the numeric data read from the excel file (without the row of column titles).

disp('File Loaded successfully!');

% Get individual tracks.

% List the points at which the trajectory number first appears:
[A,I,J] = unique(traj,'first'); % [A,I,J] = UNIQUE(traj,'first') returns the vector I to index the first occurrence of each unique value in traj.  
% A has the same values as in traj but with no repetitions. A will also be sorted.
% List the points at which the trajectory number last appears:
[A,Y,Z] = unique(traj,'last'); % UNIQUE(traj,'last'), returns the vector Y to index the last occurrence of each unique value in traj.

% Get the number of tracks (no. of different trajectories):
numtracks = numel(A);


% Create tracks structure:
tracks(1:numtracks) = struct('trajNumber',[],'xvalues0',[],'yvalues0',[],'xvalues',[],'yvalues',[],'intensity',[],'I0_IspotFit',[],'sigma_IspotFit',[],'msd_unavg',[],'frame',[],'timeabs',[],'timerel',[],'numel',[],'minNumPointsInTraj',[],'deltaTime',[],'msd',[],'errorMsd',[],'errorMsdRelPercent',[],'disp',[],'SNR',[],'bg_noise_offset_afterBGsubtract',[],'BgNoiseStd',[],'IbgAvg',[],'IinnerTot',[],'rsqFit',[]);

del = []; % initialise for later.

for i=1:numtracks 
    
    % i
    
    a = I(i); % index for starting point in trajectory.
    b = Y(i); % index for ending point in trajectory.
    
    % Delete tracks that are less than minPointsTraj data points long:
    if b-a+1 >= minPointsTraj  % Only analyse tracks which have at least "minPointsTraj" points (frames) in them (5, or 15, e.g.).
    
    data{i} = NUMERIC(a:b,:);
    % tracks(i).XLS.track_index = A(i);
    tracks(i).trajNumber = A(i);
    % all values in pixels.
    tracks(i).xvalues0 = data{i}(1:end,ID.CentreX); % original xvalues in image (used later for plotting traj on image).
    tracks(i).yvalues0 = data{i}(1:end,ID.CentreY); % original xvalues in image (used later for plotting traj on image).
    % Set origin to zero:
    tracks(i).xvalues = tracks(i).xvalues0 - (tracks(i).xvalues0(1)); % xvalues relative to the first one in the trajectory.
    tracks(i).yvalues = tracks(i).yvalues0 - (tracks(i).yvalues0(1)); % % yvalues relative to the first one in the trajectory.
    tracks(i).intensity = data{i}(1:end,ID.IspTot); % Note: we want to use IspTot (integrated spot intensity, bgnd subtracted) and not I0Fit. Isabel.    
    tracks(i).I0_IspotFit = data{i}(1:end,ID.I0Fit); % I0 from gaussian fit to intensity of spot
    tracks(i).sigma_IspotFit = data{i}(1:end,ID.SigmaFit); % sigma from gaussian fit to intensity of spot, width of spot in pixels.
    tracks(i).msd_unavg = tracks(i).xvalues.^2+tracks(i).yvalues.^2; % squared displacement from the origin: x^2 + y^2.
    tracks(i).frame = data{i}(1:end,ID.FrameNumber); % frame number.
    tracks(i).timeabs = data{i}(1:end,ID.FrameNumber).*tsamp; % tsamp is the time between frames.
    tracks(i).timerel = tracks(i).timeabs-tracks(i).timeabs(1); % Set the first frame analysed as time zero reference (not used for now). 
    tracks(i).numel = b-a+1; % Number of points in the track. Isabel: it used to be b-a, I changed it to b-a+1.
    tracks(i).minNumPointsInTraj = minPointsTraj;
    tracks(i).SNR = data{i}(1:end,ID.SNR); % Signal to noise ratio.
    tracks(i).bg_noise_offset_afterBGsubtract = data{i}(1:end,ID.bg_noise_offset_afterBGsubtract); %
    tracks(i).BgNoiseStd = data{i}(1:end,ID.BgNoiseStd); % Background noise standard deviation.
    tracks(i).IbgAvg = data{i}(1:end,ID.IbgAvg); % Intensity in background region, on average and per pixel.
    tracks(i).IinnerTot = data{i}(1:end,ID.IinnerTot); % Total intensity in circular signal mask (no bgnd subtraction).
    tracks(i).rsqFit = data{i}(1:end,ID.rsqFit); % r-square of Gaussian fit of spot intensity vs distance to spot centre for all pixels in spot signal circular mask.
    tracks(i) = getDisplacement(tracks(i),tsamp); % calculate msd and its error and add it to result structure.
    
    else
        % save indices to delete later:
        del(i) = i;     
    end
    
end

% Delete tracks which were too short: 
tracks(find(del))=[];



% Examples of how to run this function:
% Go to the directory where the .xls file which contains the trajectory data is and then for example, make a FILELIST:
% FILELIST = dir('*.xls')
% Then to analyse the first excel file in that list do:
% output = analyseTraj(FILELIST(1).name,1); )
% output is a structure with fields XLS
%     trajNumber: traj number (original one)
%     xvalues: vector with xvalues 
%     yvalues: vector with yvalues
%     intensity: vector with integrated spot intensities after background subtraction
%     msd_unavg 
%     frame : vector of frame numbers.
%     timeabs: vector of absolute times (equal to frame number only if tsamp=1.), in seconds if appropriate tsamp (time between frames in seconds) is chosen.
%     timerel: same as timeabs but all times relative to that of the first data point in trajectory. So first one is always zero.
%     numel : number of points in the track
%     deltaTime: vector of delta t values: time differences (in frames).
%     msd: vector of average mean square displacements for each delta t.
%     errorMsd: vector of absolute errors of the msd values.
%     errorMsdRelPercent: vector of relative errors of the msd values in percentage.
%     disp : cell structure: complete list of pair-wise displacements (x and y), and delta
%        time before averaging for each delta t (see getDisplacement).
%        eg: output.disp{1} is a cell structure with fields 'xdiff', 'ydiff' and
%        'tdiff' (see getDisplacement). Eg, output.disp{1}.tdiff is the
%        delta time for the differences for the first delta t,
%        output.disp{2} is the second set, etc...
%     
% Example of how to produce plots:    
% subplot(2,2,1);plot(output(2).timeabs,output(2).xvalues);title('x vs frame number');...
% subplot(2,2,3);plot(output(2).timeabs,output(2).yvalues);title('y vs frame number');...
% subplot(2,2,2);plot(output(2).timeabs,output(2).intensity,'g--x');title('Intensity vs frame number');...
% subplot(2,2,4);plot(output(2).deltaTime,output(2).msd,'r--x');title('Mean Square Displacements vs Delta t (frames)');
