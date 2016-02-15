% paramsForShowTrajAnalysis2.m

% Parameters for function ShowTrajAnalysis2.m:

%% PARAMETERS

% PARAMETER r:
% Number of frames to average over when calculating a frame average in
% order to get a Signal Mask to distinguish cell region from background
% region. Note it has to be at most equal to the total no. of frames minus
% start_frame:
% CHECK: check if only a few first frames (20) needed or more:
r = numFrames-start_frame; % number of frames to average over, starting from start_frame. (20, numFrames-start_frame)

% CHECK:
% PARAMETERs for default cell region around trajectory (for calculation of local cell coordinates):
xsize_default_cell_region = 100; % horizontal size of region around cell. (70-140) 
ysize_default_cell_region = 100; % horizontal size of region around cell. (70-140)

% PARAMETER: 
% to find out initial intensity we fit only the first
% "Npoints_for_I_line_fit" points of the Intensity vs time trace:
Npoints_for_I_line_fit = 30; % (30).
% Note that the min. number of points in the trajectory is given by
% minPointsTraj (usually 15) in showManyTrajAnalysis.m.
% If the no. of points in this traj is lower than the one above:
if analysedAllTraj(n_traj).numel < Npoints_for_I_line_fit
    Npoints_for_I_line_fit = analysedAllTraj(n_traj).numel;
end

% PARAMETERs for Chung-Kennedy filter: 'Wfilter' and 'Rfilter' are the window
% size and weighting exponent. See 'chungKennedyFilter.m'.
% Wfilter = input('Enter width of fiter window in no. of points: '); % request user input
% Rfilter = input('Enter filter weighting exponent: '); % request user input
% or:
Wfilter = 3; % (default 3-5) use for quick analysis
Rfilter = 1; % use for quick analysis
% Error control: 
% Window size for Chung-Kennedy filter cannot be larger than number of data points minus one:
Wfilter = min(Wfilter,(analysedAllTraj(n_traj).numel-1)); % choose the lowest.
    
% PARAMETER mobility_rsq_limit: value of rsquare of fit above which we
% accept that the msd versus delta-time fit is either Brownian diffusion
% (linear) or confined diffusion:
mobility_rsq_limit = 0.4;
% Maximum timeconstant in seconds from confined-trajectory fit (if time
% constant is too large, the fit is actually linear...) for trajectory to
% be labelled as confined diffusion:
max_conf_timeconst = 100; 
% Note: the guesses for the fits of the msd to a line or to a saturating
% curve are given later and their values can make the fits fair or not.

% PARAMETER nbins: number of bins in histogram of intensity pair-wise differences.
% See fullPwD.m.
nbins = 500; % (50, 200)

% PARAMETER x_limit_spectr: max Intensity step size to plot as max value in
% x-axis of power spectra for each spot...
x_limit_spectr = 3000;

% PARAMETERS for flags:
% Max number of frames in track for it to be flagged as "short" track:
max_NumFramesForShortTrack = 10; 
% Number of frames in track for it to be flagged as "long" track:
min_NumFramesForLongTrack = 50; 
% Number of frames in track for it to be flagged as "very long" track:
min_NumFramesForVeryLongTrack = 120; 
% Used for flaging trajectories which are good for extrapolating and
% obtainin the initial intensity Istart at the first frame in the track.
% Maximum no. of frames away from TimeOrigin frame of sequence:
max_framesAwayFromTimeOrigin = 10;
% For good exponential fit flags:
min_rsq_forGoodExpFit = 0.7;
max_tau_relativeError_forGoodExpFit = 30; % as percentage.

% Intensity level close to background level (we consider only ~12 steps
% of photobleaching above background). Flag tracks with low enough
% intensity levels above background to try and determine stoichiometry
% later.
% We flag up a track if it has points with intensity below these:
lowEnoughI_limit_bottom = 10000; % for bottom channel, GFP, green.
lowEnoughI_limit_top = 3000; % for top channel, mCherry, red.

