function spot_results = FindTrajects(image_label,start_frame,end_frame)
% 
% Created by Isabel Llorente-Garcia, August 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function takes an image_label such as '513', '490', etc... 
% which corresponds to a certain .sif image sequence in current folder,
% and then finds all trajectories of the bright fluorescent spots in time.
% 
% Example of how to use this function: 
% for image "Heiko_Thu Jun 24 2010_554.sif" in current folder:
% [numFrames2 frame_Ysize2 frame_Xsize2 image_data2] = extract_image_sequence_data('554'); 
%
% Reads .sif image sequence data and then for each frame:
% finds candidate spots,
% joins them to spots found and accepted on previous frame,
% eliminates coincidences in those candidates (points closer than 1 pixel),
% within all candidates, it finds actual spot centres through iterative Gaussian masking,
% and accepts only found spot centres with clipping flag equal to zero, constrained width and large enough SNR,
% eliminates coincidences from previously accepted found spot centres.
% Resulting final accepted spot centres are saved in the structure array spot_final.
% Then link spots into trajectory segments: For second frame do differently
% and check only spots in 1st frame and compare to spots in second frame.
% For rest of frames (for frame k): A) first check loose spots (TrajNumber=0) two frames
% ago (k-2) and compare to spots in current frame (k).
% B) then check all spots in previous frame (k-1) and compare to spots in
% current frame (k).
% To decide on the best asignment of pairs of spots (link), we build up
% matrices of pair-wise distances, ratio of intensities and ratio of sigmas
% of all pairs of spots in the two frames being compared, and take the
% winning asignment as that with the smallest pairwise distance.
% The intensity ratio and ratio of sigmas has to be within certain min and
% max bounds too.
%
% start_frame and end_frames are the frames through which the loop runs to
% find centres of bright spots.
%
% Output:
% The output, spot_results is a cell array with two elements:
% spot_results = {params, spot_final};
% The first element in the cell array, spot_results{1}, contains all parameters used to run
% the function: "params" (this is a structure array itself).
% The second element in the cell array, spot_results{2}, contains the track
% results: "spot_final".
% "spot_final", is itself a structure array with end_frame by L elements,
% each of which is a structure with fields: CentreX, CentreY, IspTot,
% Sigma, IbgAvg, IbgTot, IinnerTot, ClipFlag, FrameNumber and SpotNumber.
% Along the dimension L, we have the different bright spots found on each
% frame (L will usually be the number of spot centres found in frame_start,
% or the largest number of spots ever accepted).
% For example, spot_final(100,8) is a structure with the above fields
% corresponding to the eighth found spot on frame 100.
%
% Example of params: spot_results{1}:
%
%             image_label: '498'
%             start_frame: 7
%               end_frame: 500
%                       r: 5
%      max_num_candidates: 2000
%      subarray_halfwidth: 8
%     inner_circle_radius: 5
%        gauss_mask_sigma: 2
%            sigmaFit_min: 2
%            sigmaFit_max: 4
%                 SNR_min: 2
%                 rsq_min: 0.2000
%                 deflate: 1
%          d_coincid_cand: 1
%         d_coincid_found: 1
%                d_01_max: 5
%           Iratio_01_min: 0.5000
%           Iratio_01_max: 3
%       SigmaRatio_01_min: 0.5000
%       SigmaRatio_01_max: 2
%                d_02_max: 5
%           Iratio_02_min: 0.5000
%           Iratio_02_max: 3
%       SigmaRatio_02_min: 0.5000
%       SigmaRatio_02_max: 2
%                     rej: 10000
%               file_name: 'Heiko_Thu Jun 24 2010_498.sif'
%               numFrames: 500
%             frame_Ysize: 512
%             frame_Xsize: 512
%
% Example of how to call this function:
% s0 = FindTrajects('470',5,500); uses
% image sequence 470 and finds trajectories in it, for frames 5 to 500.
% Another example: s1 = FindTrajects('490',100,110);.
% -------------------------------------------------------
% Note: for reading .sif files you need 'read_sif_data_direct.m',
% 'GetAndorSifSize', etc. which are currently in path:
% 'C:\Isabel\myMatlabFiles\IO_Input\'.
%
% NOTE: before running this function you should move into the directory which
% contains the image sequence data (labelled by 'image_label').


%% DEFINITIONS and PARAMETERS:

% I have saved all parameter values in the file "paramsForFindTrajects.m"
% in the current directory. Just calling the name of that file loads the
% parameter values into the workspace:
paramsForFindTrajects
% In this way, we can save different parameter sets for different data sets
% in an easy manner and run two Matlabs at the same time working with different parameter sets.

% %
% % define data (.sif image) directory:
% % dir_data = 'Z:\Leake\Heiko Data\';
% % dir_data = 'C:\Isabel\ExperimData\HeikoData\';
% % Print data directory on command window to guide user:
% disp(' ') % empty line
% disp(['The data directory (.sif images) is: ',cd]) % display current directory.
% 
% % Save input parameters to "params":
% params.image_label = image_label;
% params.start_frame = start_frame;
% params.end_frame = end_frame;
% 
% % Number of frames to average over when calculating a frame average in
% % order to get a Signal Mask to distinguish cell region from background
% % region:
% r = 5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
% % Save parameters to results as structure "params":
% params.r = r;
% 
% % If you don't want to use the cell mask to exclude spots found outside
% it, make the following parameter equal to 1 (then spots can be accepted anywhere on image):
% doNotUseCellMaskToTellSpotsRegion = 0;
% params.doNotUseCellMaskToTellSpotsRegion = doNotUseCellMaskToTellSpotsRegion;
%
% % Maximum number of candidate spots (if eg. 260000 candidate spots are
% % found, we get an error in function pdist: "Distance matrix has more
% % elements than the maximum allowed size in MATLAB"), hence, we limit the
% % max number of candidate spots:
% max_num_candidates = 2000; % should be around 200.
% % Save parameters to results as structure "params":
% params.max_num_candidates = max_num_candidates;
% 
% % PARAMETERS for finding spot centres (see findSpotCentre1frame.m, inputs to the function):
% % The total integrated spot intensity is a bgnd corrected one, inside a
% % circular mask of radius inner_circle_radius.
% subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
% % ROI, typically a square of size 17x17 pixels. 
% inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray. 
% gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
% guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightspot intensity (default = 3).
% % Save parameters to results as structure "params":
% params.subarray_halfwidth = subarray_halfwidth;
% params.inner_circle_radius = inner_circle_radius;
% params.gauss_mask_sigma = gauss_mask_sigma;
% params.guess_sigma_Fit = guess_sigma_Fit;
% 
% % PARAMETERS for deciding if we accept a spot centre found by the function
% % findSpotCentre1frame or not:
% sigmaFit_min = -inner_circle_radius; % minimum acceptable sigma of gaussian fit to spot, in pixels (2) (-3).
% sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4) (3).
% SNR_min = 1.4; % minimum acceptable signal-to-noise ratio (at least 2) (SNR as defined in findSpotCentre1frame.m).
% rsq_min = 0.1; % minimum acceptable r-square value (0.2) (goodness of gaussian fit to spot).
% % Save parameters to results as structure "params":
% params.sigmaFit_min = sigmaFit_min;
% params.sigmaFit_max = sigmaFit_max;
% params.SNR_min = SNR_min;
% params.rsq_min = rsq_min;
% 
% % PARAMETER to decide if deflation method is applied or not (subtract each
% % found spot to allow detection of weaker spots):
% deflate = 1;
% % Save parameters to results as structure "params":
% params.deflate = deflate;
% 
% % PARAMETERS for eliminating coincident spots:
% d_coincid_cand = 1; % distance (in pixels) for eliminating coincidences in spot candidates.
% d_coincid_found = 1; % distance for eliminating coincidences in found spot centres.
% % Save parameters to results as structure "params":
% params.d_coincid_cand = d_coincid_cand; % distance for eliminating coincidences in spot candidates.
% params.d_coincid_found = d_coincid_found; % distance for eliminating coincidences in found spot centres.
% 
% % PARAMETERS for building trajectories:
% % For linking spots in current and previous frames:
% d_01_max = 5; % max distance in pixels between spot centres in current and previous frames, for linking them into a trajectory (5).
% Iratio_01_min = 0.5; % min ratio of total spot intensities (after bgnd subtraction) (0.5).
% Iratio_01_max = 3; % max ratio of total spot intensities (after bgnd subtraction) (frame k-1/frame k) (large enough value (3) to account for blinking).
% SigmaRatio_01_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit) (0.5).
% SigmaRatio_01_max = 2; % max ratio of spot width (sigma of Gaussian fit) (2).
% % Save parameters to results as structure "params":
% params.d_01_max = d_01_max; 
% params.Iratio_01_min = Iratio_01_min; 
% params.Iratio_01_max = Iratio_01_max; 
% params.SigmaRatio_01_min = SigmaRatio_01_min;
% params.SigmaRatio_01_max = SigmaRatio_01_max; 
% 
% % For linking loose spots in current frame and 2 frames ago (jump of 1 frame in trajectory):
% d_02_max = 5; % max distance in pixels between spot centres in current frame and 2 frames ago. (default: 5)
% Iratio_02_min = 0.5; % min ratio of total spot intensities (after bgnd subtraction).
% Iratio_02_max = 3; % max ratio of total spot intensities (after bgnd subtraction).
% SigmaRatio_02_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit).
% SigmaRatio_02_max = 2; % max ratio of spot width (sigma of Gaussian fit).
% % Save parameters to results as structure "params":
% params.d_02_max = d_02_max; 
% params.Iratio_02_min = Iratio_02_min; 
% params.Iratio_02_max = Iratio_02_max; 
% params.SigmaRatio_02_min = SigmaRatio_02_min; 
% params.SigmaRatio_02_max = SigmaRatio_02_max;
% 
% % Use a very large number (larger than image size in pixels) for rejected asignments:
% rej = 10000;
% params.rej = rej; % Save parameters to results as structure "params".
% 
% % Parameter to exclude a region from accepting spots (see later on, lines 322 and 507):
% exclude_region = 0;
% % Note that when exclude_region is 1, at the moment is it set for an image
% % with two channels, where spots are excluded also from around a horizontal
% % line centred on the image.
% exclude_region_width = subarray_halfwidth; % exclude edges of image, only look at spots with a centre far enough from image edges.
% params.exclude_region = exclude_region;
% params.exclude_region_width = exclude_region_width;

% Alternative way of selecting image sequence file:
% uigetfile opens a file dialog box to choose data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% strcat('data (.sif image):','  ',path_data,file_data)
% open a file dialog box to choose analysis file:
% [file_analysis,path_analysis] = uigetfile({'*.xls'}, 'Chose analysis file (trajectory):');
% strcat('analysis file (.xls trajectory):','  ',path_analysis,file_analysis)

disp(' ') % empty line
disp(['The start frame for finding bright spot trajectories will be ',num2str(start_frame)]) % start_frame is an input.
disp(['The end frame for finding bright spot trajectories will be ',num2str(end_frame)]) % end_frame is an input.
% -----------------------------------------------------


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------

% Save to parameters:
params.file_name = image_path;
params.numFrames = numFrames;
params.frame_Ysize = frame_Ysize;
params.frame_Xsize = frame_Xsize;


%% Calculate Signal Mask to distinguish cell region from background region:
% Use frame average (of first frames only) to calculate signal mask. 

% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul = zeros(frame_Ysize,frame_Xsize);

% r is the number of frames to average over, starting from start_frame.
% See PARAMETERS section.

for k = start_frame:start_frame+r  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul = frame_accumul + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg = frame_accumul/(r+1);

frame_avg_Gray = mat2gray(frame_avg); % The input to function "getCellMaskAndBoundary" needs to be a grayscale image:

% Get SignalMask to know where cells are, to distinguish cells from background:
[SignalMask CellBoundaryMask] = getCellMaskAndBoundary(frame_avg_Gray); 
% SignalMask is a matrix with 1 at positions where cells are and 0 at background.

% Or use getCellMaskAndBoundary2(frame,local_region), 
% with local_region = [xleft xright ytop ybottom].
% [SignalMask CellBoundaryMask] = getCellMaskAndBoundary2(frame_avg_Gray,[1 frame_Xsize round(frame_Ysize/2) frame_Ysize]); % Use bottom half of image only.
% Using a local region for thresholding and finding the cell mask, makes this function much faster,
% since only spots within that cell mask will be considered.

% CHECK: uncomment the following line if you don't want to use a signal
% mask. Detected spots don't need to be within a signal mask: 

if doNotUseCellMaskToTellSpotsRegion ==1
    SignalMask = ones(size(frame_avg_Gray,1),size(frame_avg_Gray,2));
end

%% Obtain candidate bright spots for start_frame (first frame), and find spot centres for those:

frame = image_data(start_frame).frame_data; % extract matrix data for first frame.
frame = double(frame);

disp(['frame number: ',num2str(start_frame)]) % print frame number to Command Window.

% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos (used in future sections):
[Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);
% Note that the image thresholding occurrs in two halves: separating top and bottom halves.
% Find candidate-bright-spots on first frame:
frame_Gray = mat2gray(frame); % The input to function "findCandidateSpots" needs to be a grayscale image:

[candidate_spotsX_00 candidate_spotsY_00] = findCandidateSpots(frame_Gray,2); % Second input: use method 2, which seems to work better.
% See C:\Isabel\myMatlabFiles\findCandidateSpots.m.
% candidate_spotsX and candidate_spotsY are two column vectors of the same
% length containing the x and y coordinates of the candidate bright spots found on the image.
% They contain integer numbers: coordinates or pixel numbers which give
% position on image plane.

% Reject candidate spots outside cell region (to speed up algorithm):
candidate_spotsX_0 = []; % initialise empty vectors before loop.
candidate_spotsY_0 = []; 
for nn = 1:length(candidate_spotsX_00)
   % Only use candidates for which there is a 1 in the SignalMask image:
   if SignalMask(candidate_spotsY_00(nn),candidate_spotsX_00(nn))==1
       candidate_spotsX_0 = [candidate_spotsX_0; candidate_spotsX_00(nn)];
       candidate_spotsY_0 = [candidate_spotsY_0; candidate_spotsY_00(nn)];
   end
end

disp(['no. of new candidate spots on start frame: ',num2str(length(candidate_spotsX_0))])

% Error control:
    % Limit the max number of candidate spots (if eg. 260000 candidate spots are
    % found, we will get an error in function pdist: "Distance matrix has more
    % elements than the maximum allowed size in MATLAB").
    % Select only the first max_num_candidates then.
    if length(candidate_spotsX_0) > max_num_candidates
        candidate_spotsX_0 = candidate_spotsX_0(1:max_num_candidates);
        candidate_spotsY_0 = candidate_spotsY_0(1:max_num_candidates);
        disp(['NOTE!! no. of candidate spots has been limited to ',num2str(max_num_candidates)])
    end

% % Check graphically:
% imshow(frame_Gray,[]);
% hold on;
% plot(candidate_spotsX_0,candidate_spotsY_0,'*');
% figure;

% Find spot centres and decide if we accept them or not:
n =1; % Initialise index n (index for accepted spot centres which have a clipping flag equal to zero):
frame_to_search = frame; % Initialise frame to search for spot centres.

% Find spot centre through iterative masking:
for m = 1:size(candidate_spotsX_0,1) % loop through all candidate spots.
    % Now find centre of bright spot using function findSpotCentre1frame:
    % use candidate spots as initial estimates and then iterate to find spot centre.
    % Image subarray ROI is a square of size 17x17 pixels (halfwidth is
    % 8 pixels), inner circular mask that moves inside the fixed 17x17
    % square has a radius of 5 pixels and the applied Gaussian
    % mask has a sigma of 2 pixels:
    spot_result = findSpotCentre1frame(frame_to_search,candidate_spotsX_0(m),candidate_spotsY_0(m),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,guess_sigma_Fit);
    spot_result.FrameNumber = start_frame; % Add new field containing frame number (time) to result structure.
    
    
    if (spot_result.ClipFlag == 0 && spot_result.noConverge == 0 && ...
            spot_result.SigmaFit <= sigmaFit_max && ...
            spot_result.SigmaFit >= sigmaFit_min && ...
            spot_result.SNR >= SNR_min &&...
            spot_result.rsqFit >= rsq_min) &&...
            (exclude_region == 0 || (exclude_region ==1 && (spot_result.CentreY <(frame_Ysize/2-exclude_region_width) || spot_result.CentreY >(frame_Ysize/2+exclude_region_width))))
            % Only accept and save result of found spot if clipping flag =0 and if values of sigmaFit, signal to noise and rsquare of fit are acceptable.
            spot_result.SpotNumber = n; % Add new field containing spot number to result structure.
            spot_final(start_frame,n) = spot_result; % store "good" found spots.
            % This is also saved in the final result spot_final, structure array.
            % first index is for frame number, second index is for spot number.
            
            %-------------------------------
            if deflate==1 % see parameter section at the beginning.
                % "Deflation" process: subtract from raw frame image the corresponding
                % Gaussian fit of each found and accepted spot before finding next spot centre (enables acceptance of dimmer spots).
                %
                % Matrices containing the x and y positions in the image frame: Xpos and Ypos
                % Xpos is a matrix of the same size as frame, containing x values for all pixels and similarly for Ypos.
                % Calculate Xpos, Ypos at the beginning: [Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);
                % Parameters of Gaussian fit of previously accepted spot:
                x_fit = spot_final(start_frame,n).CentreX;
                y_fit = spot_final(start_frame,n).CentreY;
                I_fit = spot_final(start_frame,n).I0Fit;
                sigma_fit = spot_final(start_frame,n).SigmaFit;
                % deflated frame (frame with found spot subtracted):
                deflated_frame = frame_to_search - I_fit*exp(-((Xpos-x_fit).^2+(Ypos-y_fit).^2)/(2*sigma_fit^2));
                frame_to_search = deflated_frame; % update frame to search for finding next spot-centre.
                % % Graphical check of deflated frames:
                % subplot(1,2,1); imshow(frame,[]); % frame is the original frame (always the same).
                % subplot(1,2,2); imshow(deflated_frame,[]);
            end
            %-------------------------------
            
            n = n+1; % advance index n for accepted spot centres.
    end
end

% % display the number of accepted spot-centres for this frame:
disp(['no. of accepted spot centres in first frame: ',num2str(n-1)])

% convert results of found spot-centre positions to a useful form that can
% be used as input candidate-spots on the following frame:
if (n-1) == 0 % error control: if no spots were accepted.
    found_spot_CentreX = [];
    found_spot_CentreY = [];
    % I need to create the whole spot_final structure with all its fields
    % here, just in case the number of accepted spots in the first frame is
    % zero, in order not to get error: "Subscripted assignment between
    % dissimilar structures".
    % Save empty spot (we need this, otherwise if in the last frame the no. of accepted spots is 0, there will be no result spot_final(end_frame,:) and the following functions will fail).
    spot_final(start_frame,n).CentreX = [];
    spot_final(start_frame,n).CentreY = [];
    spot_final(start_frame,n).IspTot = [];
    spot_final(start_frame,n).rsqFit = [];
    spot_final(start_frame,n).SigmaFit = [];
    spot_final(start_frame,n).I0Fit = [];
    spot_final(start_frame,n).bg_noise_offset_afterBGsubtract = [];
    spot_final(start_frame,n).BgNoiseStd = [];
    spot_final(start_frame,n).IbgAvg = [];
    spot_final(start_frame,n).IbgTot = [];
    spot_final(start_frame,n).SNR = [];
    spot_final(start_frame,n).IinnerTot = [];
    spot_final(start_frame,n).ClipFlag = [];
    spot_final(start_frame,n).noConverge = [];
    spot_final(start_frame,n).TrajNumber = [];
    spot_final(start_frame,n).FrameNumber = [];
    spot_final(start_frame,n).SpotNumber = []; 
else
    found_spot_CentreX = [spot_final(start_frame,:).CentreX]'; % column vector with found CentreX positions of all candidate spots.
    found_spot_CentreY = [spot_final(start_frame,:).CentreY]'; % column vector with found CentreY positions of all candidate spots.
end

% Check graphically:
figure;
imshow(frame_Gray,[]);
hold on;
plot(found_spot_CentreX,found_spot_CentreY,'o','Color','g','MarkerSize',10) % plot accepted spot centres in green.
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% -----------------------------------------------------------------------


%% Loop through selected frames:

tr =1; % initialise trajectory index.

for k = (start_frame+1):end_frame
    % to go through all frames do instead: for k = 1:length(sifData)
    
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    
    imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
    hold on;
    
    disp(['frame number: ',num2str(k)]) % print frame number to Command Window.
    
    %-------------------------------------
    % Find new candidate spots for this frame:
    frame_Gray = mat2gray(frame); % The input to function "findCandidateSpots" needs to be a grayscale image:
    
    [candidate_spotsX_00 candidate_spotsY_00] = findCandidateSpots(frame_Gray,2); % Second input: use method 2, which seems to work better.
    % the subindex "_00" in candidate_spotsX_00 indicates newly found spot
    % candidates for the current frame. On the other hand, found_spot_CentreX and found_spot_CentreY are the
    % accepted spot-centre positions coming from the previous frame.
    
    % Reject candidate spots outside cell region (SignalMask) (to speed up algorithm):
    candidate_spotsX_0 = []; % initialise empty vectors before loop.
    candidate_spotsY_0 = [];
    for nn = 1:length(candidate_spotsX_00)
        % Only use candidates for which there is a 1 in the SignalMask image:
        if SignalMask(candidate_spotsY_00(nn),candidate_spotsX_00(nn))==1
            candidate_spotsX_0 = [candidate_spotsX_0; candidate_spotsX_00(nn)];
            candidate_spotsY_0 = [candidate_spotsY_0; candidate_spotsY_00(nn)];
        end
    end

    disp(['no. of new candidate spots: ',num2str(length(candidate_spotsX_0))])
    %------------------------------------- 
    
    % Join accepted spot-centre positions from previous frame with
    % candidate spots for this frame to use them as new candidates for this frame:
    candidate_spotsX = [found_spot_CentreX; candidate_spotsX_0];
    candidate_spotsY = [found_spot_CentreY; candidate_spotsY_0];
    
%          % Plot new candidate spots in yellow and found spot centres from
%          % previous frame in cyan
%           plot(candidate_spotsX_0,candidate_spotsY_0,'+','Color','y','MarkerSize',3);
%           pause(0.5);
%           plot(found_spot_CentreX,found_spot_CentreY,'+','Color','c','MarkerSize',3);
%           pause(0.5);

    disp(['no. of initial total candidate spots: ',num2str(length(candidate_spotsX))])
    
    % Error control:
    % Limit the max number of candidate spots (if eg. 260000 candidate spots are
    % found, we get an error in function pdist: "Distance matrix has more
    % elements than the maximum allowed size in MATLAB").
    % Select only the first max_num_candidates then.
    if length(candidate_spotsX) > max_num_candidates
        candidate_spotsX = candidate_spotsX(1:max_num_candidates);
        candidate_spotsY = candidate_spotsY(1:max_num_candidates);
        disp(['NOTE!! no. of initial total candidate spots has been limited to ',num2str(max_num_candidates)])
    end
    
    %----------------------------------------------
    % Eliminate coincidences in spot candidates:
    
    [candidate_spotsX candidate_spotsY pos] = eliminateCoincidentSpots(candidate_spotsX,candidate_spotsY,d_coincid_cand);
    % see C:\Isabel\myMatlabFiles\eliminateCoincidentSpots.m
    % The function checks the distances between all pairs of points with x
    % and y coordinates candidate_spotsX and candidate_spotsY respectively,
    % and removes those points (x,y) which are closer than one pixel (distance<1) to
    % another point in the list.
    
    %     % for debugging:
    %     x2 = candidate_spotsX;
    %     y2 = candidate_spotsY;
    
    disp(['no. of total candidate spots after eliminating coincidences: ',num2str(length(candidate_spotsX))])
    
    %     % Plot all candidate spots in magenta after removing coincidences:
    %      plot(candidate_spotsX,candidate_spotsY,'+','Color','m','MarkerSize',3);
    %      pause(0.5);
    %------------------------------------------------
    
    
    n =1; % Initialise index n (index for accepted spot centres which have a clipping flag equal to zero):
    frame_to_search = frame; % Initialise frame to search for spot centres.
    
    for m = 1:size(candidate_spotsX,1) % for each frame, loop throuh all the candidate spots.
        % Now find centre of bright spot using function findSpotCentre1frame:
        % use candidate spots as initial estimates and then iterate to find spot centre.
        % Image subarray ROI is a square of size 17x17 pixels (halfwidth is
        % 8 pixels), inner circular mask that moves inside the fixed 17x17
        % square has a radius of 5 pixels and the applied Gaussian
        % mask has a sigma of 2 pixels:
        spot_result = findSpotCentre1frame(frame_to_search,candidate_spotsX(m),candidate_spotsY(m),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,guess_sigma_Fit);
        % index k is for frame number, index m is for spot number
        spot_result.FrameNumber = k; % Add new field containing frame number (time) to result structure.
        
        % accepted spot centres:
        if (spot_result.ClipFlag == 0 && spot_result.noConverge == 0 && ...
                spot_result.SigmaFit <= sigmaFit_max && ...
                spot_result.SigmaFit >= sigmaFit_min && ...
                spot_result.SNR >= SNR_min &&...
                spot_result.rsqFit >= rsq_min)&&...
                (exclude_region == 0 || (exclude_region ==1 && (spot_result.CentreY <(frame_Ysize/2-exclude_region_width) || spot_result.CentreY >(frame_Ysize/2+exclude_region_width))))
            % Only accept and save result of found spot if clipping flag =0 and if values of sigmaFit, signal to noise and rsquare of fit are acceptable.
            spot(k,n) = spot_result; % store accepted found spots in this preliminary result.
            % first index is for frame number, second index is for spot number.
            %           %--------------------------------------------------
            %           plot(spot(k,n).CentreX,spot(k,n).CentreY,'o','Color','r','MarkerSize',10); % Plot found centre spots in red
            %           %--------------------------------------------------
            if deflate==1 % see parameter section at the beginning.
                % "Deflation" process: subtract from raw frame image the corresponding
                % Gaussian fit of each found and accepted spot before finding
                % next spot centre (enables acceptance of dimmer spots).
                % Xpos and Ypos are matrices of x and y positions on image frame.
                % Parameters of Gaussian fit of previously accepted spot:
                x_fit = spot(k,n).CentreX;
                y_fit = spot(k,n).CentreY;
                I_fit = spot(k,n).I0Fit;
                sigma_fit = spot(k,n).SigmaFit;
                % deflated frame (frame with found spot subtracted):
                deflated_frame = frame_to_search - I_fit*exp(-((Xpos-x_fit).^2+(Ypos-y_fit).^2)/(2*sigma_fit^2));
                frame_to_search = deflated_frame; % update frame for finding next spot-centre.
                % % Graphical check of deflated frames:
                % subplot(1,2,1); imshow(frame,[]); % frame is the original frame (always the same).
                % subplot(1,2,2); imshow(deflated_frame,[]);
            end
            %-------------------------------
            
            n = n+1; % and advance index n for accepted spot centres.
        end
        
    end
    
    % display the number of accepted spot-centres for each frame:
    disp(['no. of accepted spot centres: ',num2str(n-1)])
    
    % % The following two lines are used together with the previous two
    % "plot" and "imshow" (commented off) lines:
    pause(0.1); % this pause is needed to give time for the plot to appear
    %    hold off;
    
    % convert results of found spot-centre positions to a useful form that can
    % be used as input candidate-spots on the following frame:
    if (n-1) == 0 % error control: if no spots were accepted.
        found_spot_CentreX = [];
        found_spot_CentreY = [];
        spot_final(k,n).SpotNumber = []; % Save empty spot (we need this, otherwise if in the last frame the no. of accepted spots is 0, there will be no result spot_final(end_frame,:) and the following functions will fail).
    else
        found_spot_CentreX = [spot(k,:).CentreX]'; % column vector with found CentreX positions of all candidate spots.
        found_spot_CentreY = [spot(k,:).CentreY]'; % column vector with found CentreY positions of all candidate spots.
        
        %-------------------------------
        % Eliminate coincidences in result of last found spots for a given frame (for distance <1):
        [found_spot_CentreX found_spot_CentreY pos_final] = eliminateCoincidentSpots(found_spot_CentreX,found_spot_CentreY,d_coincid_found);
        % see C:\Isabel\myMatlabFiles\eliminateCoincidentSpots.m
        % pos_final contains positions of selected, kept spot centres.
        %-------------------------------
        
        % Save final spots to variable final_spots:
        n=1; % index for final kept spot.
        for ii = 1:length(pos_final)
            mientras = spot(k,pos_final(ii)); % intermediate result.
            mientras.SpotNumber = n; % Add new field containing spot number to result structure.
            spot_final(k,n)=mientras; % final result structure of accepted spot centres.
            n = n+1;
        end
        % Plot found spot centres:
        pause(0.5);
        plot(found_spot_CentreX,found_spot_CentreY,'o','Color','g','MarkerSize',10) % plot final accepted spot centres in green.
        pause(0.1); % this pause is needed to give time for the plot to appear
        hold off;
        
        disp(['no. of final found spot centres after eliminating coincidences: ',num2str(length(found_spot_CentreX))])
    end
    
    
    
    %--------------------
    % LINKING SPOTS INTO TRAJECTORY SEGMENTS:
    
    % Link found and accepted spots into trajectory segments:
    
    % Trajectory index tr is initialised to 1 outside the loop through frames (k loop).
    
    % Do differently FOR SECOND FRAME (k == start_frame+1): compare only accepted spots in
    % previous and current frames:
    if k == start_frame+1 && ... % If second frame and
            (n-1)~=0 && ... % if the number of accepted spot centres is not zero and
            isempty([spot_final(k-1,:).SpotNumber])==0 && ... % at least 1 accepted spot in previous frame and
            isempty([spot_final(k,:).SpotNumber])==0 % at least 1 accepted spot in current frame.
        % There are no trajectories jet, so compare accepted spots in previous and current frames:
        N0 = max(cat(1,spot_final(k-1,:).SpotNumber));  % no. of accepted spots in previous frame.
        % Note: cat(1,spot_final(k-1,:).SpotNumber) gives a column vector with the values of SpotNumber for all non-empty accepted spots in frame k-1.
        N1 = max(cat(1,spot_final(k,:).SpotNumber));  % no. of accepted spots in current frame.
        
        % Create cell arrays with empty elements to pre-asign sizes:
        d01 = cell(N0,N1); % Note: d01 is a cell array (matrix) but d_01 below is a scalar.
        Iratio01 = cell(N0,N1); % Note: Iratio01 is a cell array (matrix) but Iratio_01 below is a scalar.
        SigmaRatio01 = cell(N0,N1); % Note: SigmaRatio01 is a cell array (matrix) but SigmaRatio_01 below is a scalar.
        
        for q0 = 1:N0 % loop though accepted spots in previous frame.
            for q1 = 1:N1 % loop though accepted spots in current frame.
                % d_01: distance between spot centres in previous and current frames:
                d_01 = sqrt((spot_final(k-1,q0).CentreX-spot_final(k,q1).CentreX)^2+(spot_final(k-1,q0).CentreY-spot_final(k,q1).CentreY)^2);
                % Iratio_01: ratio of intensities of spot centre in previous and current frames:
                Iratio_01 = spot_final(k-1,q0).IspTot/spot_final(k,q1).IspTot;
                % SigmaRatio_01: ratio of widths of spots (Gaussian fits) in previous and current frames:
                SigmaRatio_01 = spot_final(k-1,q0).SigmaFit/spot_final(k,q1).SigmaFit;
                
                %                 d_01
                %                 Iratio_01
                %                 SigmaRatio_01
                
                % Accept and save trajectory if spots in previous and
                % current frames fulfill the following conditions:
                if d_01 < d_01_max && ...  % see PARAMETERS at start of this function.
                        Iratio_01_min <= Iratio_01 && Iratio_01 <= Iratio_01_max && ...
                        SigmaRatio_01_min <= SigmaRatio_01 && SigmaRatio_01 <= SigmaRatio_01_max
                    % Asign accepted values to cell array elements to store them:
                    d01{q0,q1} = d_01; % use {} for cell arrays.
                    Iratio01{q0,q1} = Iratio_01;
                    SigmaRatio01{q0,q1} = SigmaRatio_01;
                else % rejected asignments:
                    d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                    Iratio01{q0,q1} = rej; % Use rej for asignments not accepted.
                    SigmaRatio01{q0,q1} = rej; % Use rej for rejected asignments.
                end
            end
            
%                         d01
%                         [d01{q0,:}]
            % Note that [d01{q0,:}] gives only non-empty elements of row q0
            % in the cell array d01 as a row vector, that's why we had to
            % give a numeric value rej to non-accepted asignments.
            
            % Note that if all asignments in previous step are rejected,
            % [d01{q0,:}] will be a list of rej values, and its minimum will
            % be rej.
            % If list of "linkable" spots, [d01{q0,:}], has no accepted
            % asignments (all values are rej):
            if min([d01{q0,:}]) == rej
                % Asign trajectory number 0 to the spot in the previous frame only:
                spot_final(k-1,q0).TrajNumber = 0;
            else % if there is at least one accepted asignment for a given spot in the previous frame:
                
                % Decide of all possible accepted spots (in current frame)
                % that could be linked to spot q0 in previous frame, which one is the best:
                % We take the best as the closest one to spot q0:
                q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                
%                             q1_chosen
                
                % Check if there is a better competing asignment for a given spot q1 in the current
                % frame from another spot in the previous frame.
                % Hence, check also column-wise in matrix d01 to avoid asigning a traj
                % number to a spot q1 in the current frame that had already
                % had a traj number asigned to it linking it to a different spot q0 in
                % the previous frame, which might be at a shorter distance
                % from it than the current one.
                
%                             [d01{:,q1_chosen}] % chosen column of d01 matrix of distances.
                
                % Asign trajectory numbers to structure spot_final:
                % If the found distance in that column is not the minimum one:
                if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                    spot_final(k-1,q0).TrajNumber = 0; % asign trajectory number 0 to spot in previous frame.
                else
                    spot_final(k-1,q0).TrajNumber = tr; % asign trajectory number to spot in previous frame, to spot_final structure.
                    spot_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to spot in current frame.
                    tr = tr+1; % advance trajectory-number index.
                end
            end
        end
        
        
        
    else % for FRAMES k >= start_frame+2, from third chosen frame on:
        
        % A) Compare loose spots (TrajNumber is 0) two frames ago (k-2)
        % to found spots in current frame (TrajNumber is []):
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % DO maybe!!
        if k ~= start_frame+1 && (n-1)~=0 && ... % If the number of accepted spot centres is not zero and
                isempty([spot_final(k-2,:).SpotNumber])==0 && ... % at least 1 accepted spot 2 frames ago and
                isempty([spot_final(k,:).SpotNumber])==0 % at least 1 accepted spot in current frame.
            
            N0 = max(cat(1,spot_final(k-2,:).SpotNumber));  % no. of accepted spots 2 frames ago.
            N1 = max(cat(1,spot_final(k,:).SpotNumber));  % no. of accepted spots in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d02 = cell(N0,N1); % Note: d02 is a cell array (matrix) but d_02 below is a scalar.
            Iratio02 = cell(N0,N1); % Note: Iratio02 is a cell array (matrix) but Iratio_02 below is a scalar.
            SigmaRatio02 = cell(N0,N1); % Note: SigmaRatio02 is a cell array (matrix) but SigmaRatio_02 below is a scalar.
            
            for q0 = 1:N0 % loop though loose accepted spots 2 frames ago.
                if spot_final(k-2,q0).TrajNumber == 0 % only for loose (unlinked) spots (TrajNumber=0) two frames ago (so only rows q0 in matrix d01 which have unlinked spots will fill up).
                    for q1 = 1:N1 % loop though accepted spots in current frame.
                        % d_01: distance between spot centres in previous and current frames:
                        d_02 = sqrt((spot_final(k-2,q0).CentreX-spot_final(k,q1).CentreX)^2+(spot_final(k-2,q0).CentreY-spot_final(k,q1).CentreY)^2);
                        % Iratio_01: ratio of intensities of spot centre in previous and current frames:
                        Iratio_02 = spot_final(k-2,q0).IspTot/spot_final(k,q1).IspTot;
                        % SigmaRatio_01: ratio of widths of spots (Gaussian fits) in previous and current frames:
                        SigmaRatio_02 = spot_final(k-2,q0).SigmaFit/spot_final(k,q1).SigmaFit;
                        
                        %                         d_02
                        %                         Iratio_02
                        %                         SigmaRatio_02
                        
                        % Accept and save trajectory if spots in previous and
                        % current frames fulfill the following conditions:
                        if d_02 < d_02_max && ...  % see PARAMETERS at start of this function.
                                Iratio_02_min <= Iratio_02 && Iratio_02 <= Iratio_02_max && ...
                                SigmaRatio_02_min <= SigmaRatio_02 && SigmaRatio_02 <= SigmaRatio_02_max
                            % Asign accepted values to cell array elements to store them:
                            d02{q0,q1} = d_02; % use {} for cell arrays.
                            Iratio02{q0,q1} = Iratio_02;
                            SigmaRatio02{q0,q1} = SigmaRatio_02;
                        else % rejected asignments:
                            d02{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                            Iratio02{q0,q1} = rej; % Use rej for asignments not accepted.
                            SigmaRatio02{q0,q1} = rej; % Use rej for rejected asignments.
                        end
                    end
                    
%                                         d02
%                                         [d02{q0,:}]
                    % Note that [d02{q0,:}] gives only non-empty elements of the cell array d01 as a row vector.
                    
                    % Note that if all asignments in previous step are rejected,
                    % [d02{q0,:}] will be a list of rej values, and its minimum will be rej.
                    % If list of "linkable" spots, [d02{q0,:}], has no accepted asignments (all values are rej):
                    if min([d02{q0,:}]) == rej
                        % Asign trajectory number 0 to the spot two frames ago only:
                        spot_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        
                    else % if there is at least one accepted asignment for a given spot two frames ago:
                        
                        % Decide of all possible accepted/saved spots (in current frame)
                        % that could be linked to spot q0 in previous frame (of all possible asignments), which one is the best:
                        % We take the best as the closest one to spot q0:
                        q1_chosen = find([d02{q0,:}] == min([d02{q0,:}])); % find position of the minimum pair-wise distance.
                        
%                                             q1_chosen
                        
                        % Check if there is a better competing asignment for a given spot q1 in the current
                        % frame from another spot q0 two frames ago.
                        % Hence, check also column-wise in d01 to avoid asigning a traj
                        % number to a spot q1 in the current frame that had already
                        % had a traj number asigned to it linking it to a different spot q0 two frames ago which might be at a shorter distance
                        % from it than the current one.
                        
%                                             [d02{:,q1_chosen}] % chosen column of d01 matrix of distances.
                        
                        % Asign trajectory numbers to structure spot_final:
                        % If the found distance in that column is not the minimum one:
                        if q0 ~= find([d02{:,q1_chosen}] == min([d02{:,q1_chosen}]));
                            % Asign trajectory number 0 to the spot two frames ago only:
                            spot_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        else
                            spot_final(k-2,q0).TrajNumber = tr; % asign trajectory number to spot two frames ago, to spot_final structure.
                            spot_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to spot in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
        end
        
        
        % B) Compare loose spots (TrajNumber is []) and trajectories (TrajNumber is >0)
        % in previous frame (k-1) to found spots in current frame:
        if (n-1)~=0 && ... % If the number of accepted spot centres is not zero and
                isempty([spot_final(k-1,:).SpotNumber])==0 && ... % at least 1 accepted spot in previous frame.
                isempty([spot_final(k,:).SpotNumber])==0 % at least 1 accepted spot in current frame.
            
            % zzzzzzzzzzzzzzzzzzzzzzz
            N0 = max(cat(1,spot_final(k-1,:).SpotNumber));  % no. of accepted spots in previous frame.
            N1 = max(cat(1,spot_final(k,:).SpotNumber));  % no. of accepted spots in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d01 = cell(N0,N1);
            Iratio01 = cell(N0,N1);
            SigmaRatio01 = cell(N0,N1);
            
            for q0 = 1:N0 % loop though accepted spots in previous frame.
                for q1 = 1:N1 % loop though accepted spots in current frame.
                    % d_01: distance between spot centres in previous and current frames:
                    d_01 = sqrt((spot_final(k-1,q0).CentreX-spot_final(k,q1).CentreX)^2+(spot_final(k-1,q0).CentreY-spot_final(k,q1).CentreY)^2);
                    % Iratio_01: ratio of intensities of spot centre in previous and current frames:
                    Iratio_01 = spot_final(k-1,q0).IspTot/spot_final(k,q1).IspTot;
                    % SigmaRatio_01: ratio of widths of spots (Gaussian fits) in previous and current frames:
                    SigmaRatio_01 = spot_final(k-1,q0).SigmaFit/spot_final(k,q1).SigmaFit;
                    
                    %                     d_01
                    %                     Iratio_01
                    %                     SigmaRatio_01
                    
                    % Accept and save trajectory if spots in previous and
                    % current frames fulfill the following conditions:
                    if d_01 < d_01_max && ...  % see PARAMETERS at start of this function.
                            Iratio_01_min <= Iratio_01 && Iratio_01 <= Iratio_01_max && ...
                            SigmaRatio_01_min <= SigmaRatio_01 && SigmaRatio_01 <= SigmaRatio_01_max
                        % Asign accepted values to cell array elements to store them:
                        d01{q0,q1} = d_01; % use {} for cell arrays.
                        Iratio01{q0,q1} = Iratio_01;
                        SigmaRatio01{q0,q1} = SigmaRatio_01;
                    else % rejected asignments:
                        d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                        Iratio01{q0,q1} = rej; % Use rej for asignments not accepted.
                        SigmaRatio01{q0,q1} = rej; % Use rej for rejected asignments.
                    end
                end
                
%                                 d01
%                                 [d01{q0,:}] % last row of d01 matrix of distances.
                % Note that [d01{q0,:}] gives only non-empty elements
                % of row q0 in the cell array d01, as a row vector.
                
                % Note that if all asignments in previous step are rejected,
                % [d01{q0,:}] will be a list of rej values, and its minimum will
                % be rej.
                % If list of "linkable" spots, [d01{q0,:}], has no accepted
                % asignments (all values are rej):
                if min([d01{q0,:}]) == rej
                    
                    if isempty(spot_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                        % Asign trajectory number 0 to the spot in the previous frame only:
                        spot_final(k-1,q0).TrajNumber = 0;                       
                    end
                    
                else % if there is at least one accepted asignment for a given spot two frames ago:
                    
                    % Decide of all possible accepted/saved spots (in current frame)
                    % that could be linked to spot q0 in previous frame, which one is the best:
                    % We take the best as the closest one to spot q0:
                    q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                    
%                                     q1_chosen
                    
                    % Check if there is a better competing asignment for a given spot q1 in the current
                    % frame from another spot q0 in the previous frame.
                    % Hence, check also column-wise in d01 to avoid asigning a traj
                    % number to a spot q1 in the current frame that had already
                    % had a traj number asigned to it linking it to a different spot q0 in
                    % the previous frame which might be at a shorter distance
                    % from it than the current one.
                    
%                                     [d01{:,q1_chosen}]  % chosen column of d01 matrix of distances.
                    
                    % If the found distance in that column is not the minimum one:
                    if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                        if isempty(spot_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                            % Asign trajectory number 0 to the spot in the previous frame only:
                            spot_final(k-1,q0).TrajNumber = 0;
                        end
                    else
                        % Asign trajectory numbers to structure spot_final:
                        if spot_final(k-1,q0).TrajNumber > 0 % if point in previous frame was already part of a trajectory:
                            spot_final(k,q1_chosen).TrajNumber = spot_final(k-1,q0).TrajNumber; % asign that trajectory number to spot in current frame.
                        else % if point in previous frame was not part of a trajectory:
                            spot_final(k-1,q0).TrajNumber = tr; % asign new trajectory number to spot in previous frame.
                            spot_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to spot in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
            % zzzzzzzzzzzzzzzzzzzzzzz
        end
    end
    %-----------------
    
    
end  % loop through selected frames



%% OUTPUT OF SPOT-FINDING PROCESS: final output spot_results:
%
spot_results = {params, spot_final};
% params is a structure array containing all parameters used to run the
% function.
% spot_final, is a structure array with end_frame x L elements,
% each of which is a structure with fields:
%     'CentreX'
%     'CentreY'
%     'IspTot'
%     'rsqFit'
%     'SigmaFit'
%     'I0Fit'
%     'BgNoiseStd'
%     'IbgAvg'
%     'IbgTot'
%     'SNR'
%     'IinnerTot'
%     'ClipFlag'
%     'TrajNumber'
%     'FrameNumber'
%     'SpotNumber'
%
% Along the dimension L, we have the different bright spots found on each
% frame (L will often be the number of spot centres found in frame_start,
% it is always the largest number of spots ever accepted on one frame).
% For example, spot_final(100,8) is a structure with the above fields
% corresponding to the eighth found spot on frame 100.
%
% Note that even if we only analyse from start_frame to end_frame,
% spot_final is a list containing empty structure arrays from index 1 to
% index start_frame, and then the found spots for the analysed frames start_frame to end_frame.
%
% The result is padded to a fixed number of spot structures for reach
% frame (the maximum no. of accepted found spots of all frames), so that for a given
% frame in which less found spots have been accepted, the remaining
% elements are padded with empty structures with empty fields [].
% To check if a given spot is empty: isempty(spot_final(101,10).CentreX)
% gives 1 if field "CentreX" of tenth spot found and accepted in frame 101
% is empty (equal to []).
%
% e.g. cat(1,spot_final(100,:).SigmaFit) gives a vector column with all the
% non-empty SigmaFit values for all spot centres accepted in frame 100.
