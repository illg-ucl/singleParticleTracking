%% DEFINITIONS and PARAMETERS:

% paramsForFindTrajects.m

%
% define data (.sif image) directory:
% dir_data = 'Z:\Leake\Heiko Data\';
% dir_data = 'C:\Isabel\ExperimData\HeikoData\';
% Print data directory on command window to guide user:
disp(' ') % empty line
disp(['The data directory (.sif images) is: ',cd]) % display current directory.

% Save input parameters to "params":
params.image_label = image_label;
params.start_frame = start_frame;
params.end_frame = end_frame;

% Number of frames to average over when calculating a frame average in
% order to get a Signal Mask to distinguish cell region from background
% region:
r = 5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
% Save parameters to results as structure "params":
params.r = r;

% If you don't want to use the cell mask to exclude spots found outside it,
% make the following parameter equal to 1 (then spots can be accepted
% anywhere on image):
doNotUseCellMaskToTellSpotsRegion = 0;
params.doNotUseCellMaskToTellSpotsRegion = doNotUseCellMaskToTellSpotsRegion;

% Maximum number of candidate spots (if eg. 260000 candidate spots are
% found, we get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB"), hence, we limit the
% max number of candidate spots:
max_num_candidates = 2000; % should be around 200.
% Save parameters to results as structure "params":
params.max_num_candidates = max_num_candidates;

% PARAMETERS for finding spot centres (see findSpotCentre1frame.m, inputs to the function):
% The total integrated spot intensity is a bgnd corrected one, inside a
% circular mask of radius inner_circle_radius.
subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels. 
inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray. 
gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightspot intensity (default = 3).
% Save parameters to results as structure "params":
params.subarray_halfwidth = subarray_halfwidth;
params.inner_circle_radius = inner_circle_radius;
params.gauss_mask_sigma = gauss_mask_sigma;
params.guess_sigma_Fit = guess_sigma_Fit;

% PARAMETERS for deciding if we accept a spot centre found by the function
% findSpotCentre1frame or not:
sigmaFit_min = -inner_circle_radius; % minimum acceptable sigma of gaussian fit to spot, in pixels (2) (-3).
sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4) (3).
SNR_min = 2; % minimum acceptable signal-to-noise ratio (at least 2) (SNR as defined in findSpotCentre1frame.m).
rsq_min = 0.2; % minimum acceptable r-square value (0.2) (goodness of gaussian fit to spot).
% Save parameters to results as structure "params":
params.sigmaFit_min = sigmaFit_min;
params.sigmaFit_max = sigmaFit_max;
params.SNR_min = SNR_min;
params.rsq_min = rsq_min;

% PARAMETER to decide if deflation method is applied or not (subtract each
% found spot to allow detection of weaker spots):
deflate = 1;
% Save parameters to results as structure "params":
params.deflate = deflate;

% PARAMETERS for eliminating coincident spots:
d_coincid_cand = 1; % distance (in pixels) for eliminating coincidences in spot candidates.
d_coincid_found = 1; % distance for eliminating coincidences in found spot centres.
% Save parameters to results as structure "params":
params.d_coincid_cand = d_coincid_cand; % distance for eliminating coincidences in spot candidates.
params.d_coincid_found = d_coincid_found; % distance for eliminating coincidences in found spot centres.

% PARAMETERS for building trajectories:
% For linking spots in current and previous frames:
d_01_max = 5; % max distance in pixels between spot centres in current and previous frames, for linking them into a trajectory (5).
Iratio_01_min = 0.5; % min ratio of total spot intensities (after bgnd subtraction) (0.5).
Iratio_01_max = 3; % max ratio of total spot intensities (after bgnd subtraction) (frame k-1/frame k) (large enough value (3) to account for blinking).
SigmaRatio_01_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit) (0.5).
SigmaRatio_01_max = 2; % max ratio of spot width (sigma of Gaussian fit) (2).
% Save parameters to results as structure "params":
params.d_01_max = d_01_max; 
params.Iratio_01_min = Iratio_01_min; 
params.Iratio_01_max = Iratio_01_max; 
params.SigmaRatio_01_min = SigmaRatio_01_min;
params.SigmaRatio_01_max = SigmaRatio_01_max; 

% For linking loose spots in current frame and 2 frames ago (jump of 1 frame in trajectory):
d_02_max = 5; % max distance in pixels between spot centres in current frame and 2 frames ago. (default: 5)
Iratio_02_min = 0.2; % min ratio of total spot intensities (after bgnd subtraction).
Iratio_02_max = 3; % max ratio of total spot intensities (after bgnd subtraction).
SigmaRatio_02_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit).
SigmaRatio_02_max = 2; % max ratio of spot width (sigma of Gaussian fit).
% Save parameters to results as structure "params":
params.d_02_max = d_02_max; 
params.Iratio_02_min = Iratio_02_min; 
params.Iratio_02_max = Iratio_02_max; 
params.SigmaRatio_02_min = SigmaRatio_02_min; 
params.SigmaRatio_02_max = SigmaRatio_02_max;

% Use a very large number (larger than image size in pixels) for rejected asignments:
rej = 10000;
params.rej = rej; % Save parameters to results as structure "params".

% Parameter to exclude a region from accepting spots (see later on, lines 322 and 507):
exclude_region = 0;
% Note that when exclude_region is 1, at the moment it is set for an image
% with two channels, where spots are excluded also from around a horizontal
% line centred on the image.
exclude_region_width = subarray_halfwidth; % exclude edges of image, only look at spots with a centre far enough from image edges.
params.exclude_region = exclude_region;
params.exclude_region_width = exclude_region_width;
