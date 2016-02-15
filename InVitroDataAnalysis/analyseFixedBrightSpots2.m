function result_final = analyseFixedBrightSpots2(image_label,start_frame,end_frame,folder_label)
%
% March 2012. Created by Isabel Llorente Garcia. 
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function has a lot in common with findTrajects.m
% Same as analyseFixedBrightSpots.m but plots more results and analyses
% each spot in two different ways. Takes a lot longer to run too!!
%
% This function is used to analyse fluorescent labels which are fixed
% (glued) onto a glass surface when imaged. This means tracking is not necessary.
% I first get a frame average and find the position of the fixed bright
% spots and then monitor the fluorescence in those spots over time.
% It should be possible to extract on-off rates of blinking and histograms
% of intensity levels.
% The aim is to extract the typical intensity of a single-molecule of a
% given fluorophore.
%
% Takes ~1 minute to run for ~100 frames and around 70 spots.
%
% Inputs:
% image_label: string that labels a given sequence: eg. '470'.
% start_frame: number for first frame to look at: eg. 5.
% end_frame: number for last frame to look at: eg. 100.
% User input required for ROI (region of interest) limits: eg. enter  [1 500 250 500]
% for [xleft xright ytop ybottom] ROI limits.
% folder_label: label of choice to be added to the end of the output folder name.
%
% Output: "result_final" is a cell array with two elements:
% result_final = {results params};
% To obtain "results" do: result_final{1}.
% To obtain "params" do: result_final{2}.
%
% "results" is a structure array (containing data for found spots), each with fields: 
%     spot_number
%     Xcentre
%     Ycentre
%     frame_number
%     integrated_Ispot
%     integrated_Ispot_bgndSubtracted
%     integrated_Ispot_corrected
%     filtered_Ispot
%     binCentres_Ipwd
%     freqCounts_Ipwd
%     power_spectrum_X
%     power_spectrum_Y
%     power_spectrum_peaks_X
%     power_spectrum_peaks_Y
%     power_spectrum_X_neg
%     power_spectrum_Y_neg
%     power_spectrum_peaks_X_neg
%     power_spectrum_peaks_Y_neg
%
% "params": is a structure containing as fields the different parameters
% used for the analysis:
%     'image_label'
%     'folder_label'
%     'start_frame'
%     'end_frame'
%     'max_num_candidates'
%     'sigmaFit_min'
%     'sigmaFit_max'
%     'SNR_min'
%     'rsq_min'
%     'deflate'
%     'square_halfsize'
%     'Wfilter'
%     'Rfilter'
%     'nbins'
%     'x_limit_spectr'
%     'image_path'
%     'numFrames'
%     'frame_Ysize'
%     'frame_Xsize'
%     'xleft_ROI'
%     'xright_ROI'
%     'ytop_ROI'
%     'ybottom_ROI'
%     'bgnd_spot_x'
%     'bgnd_spot_y'
%     'Ibgnd' - vector
%     'Ibgnd_filtered' - vector
%     'Ibgnd_filtered_normalised' - vector
    
%% DEFINITIONS and PARAMETERS:
%
% Print data directory on command window to guide user:
disp(' '); % empty line
disp(['The data directory is: ',cd]) % display current directory.

% Save inputs as parameters for output:
params.image_label = image_label;
params.folder_label = folder_label;
params.start_frame = start_frame;
params.end_frame = end_frame;

% Maximum number of candidate spots (if eg. 260000 candidate spots are
% found, we get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB"), hence, we limit the
% max number of candidate spots:
params.max_num_candidates = 10000; % (2000)

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
params.sigmaFit_min = -3.5; % minimum acceptable sigma of gaussian fit to spot, in pixels (2).
params.sigmaFit_max = 3.5; % maximum acceptable sigma of gaussian fit to spot, in pixels (4).
params.SNR_min = 2; % minimum acceptable signal-to-noise ratio (at least 2) (SNR as defined in findSpotCentre1frame.m).
params.rsq_min = 0.2; % minimum acceptable r-square value (0.2) (goodness of gaussian fit to spot).

% PARAMETER deflate: to decide if deflation method is applied or not (subtract each
% found spot to allow detection of weaker spots):
params.deflate = 1; % 1 for yes, 0 for no.

% PARAMETER square_halfsize: size of circular mask to measure intensity around position of
% each spot found and accepted (default is 5 pixels):
params.square_halfsize = 6; % (4 for Rodrigo's data, 6-8 for Oxphos data)

% PARAMETERS for Chung-Kennedy filter: 'Wfilter' and 'Rfilter' are the window
% size and weighting exponent. See 'chungKennedyFilter.m'.

% Wfilter = input('Enter width of fiter window in no. of points: '); % request user input
% Rfilter = input('Enter filter weighting exponent: '); % request user input
% or:
params.Wfilter = 3; % (default 4-5) use for quick analysis
params.Rfilter = 1; % use for quick analysis (1).

% PARAMETER points_per_bin: to decide size of bins in histogram of intensity pair-wise differences.
% See fullPwD.m.
params.nbins = 100; % (100)

% PARAMETER x_limit_spectr: max Intensity step size to plot as max value in
% x-axis of power spectra for each spot...
params.x_limit_spectr = 1200; % (5000 for Rodrigo's data, 1200-3000 for Oxphos data)
%----------------------------------------

%% Create new folder for saving output figures and results:

% Make new folder (new directory) to save png figures and results for each image sequence:
result_folder_name = strcat('resultsFixedSpots',image_label,folder_label);
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(result_folder_name); % make new directory.


%% Read in the Image-Sequence data:
% (The next few sections are copied from frameAverage2.m)

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.

disp(' ') % empty line
disp(['The start frame is ',num2str(start_frame)]) % start_frame is an input.
disp(['The end frame is ',num2str(end_frame)]) % end_frame is an input.

% Save as parameters for output:
params.image_path = image_path;
params.numFrames = numFrames;
params.frame_Ysize = frame_Ysize;
params.frame_Xsize = frame_Xsize;
% --------------------------------

%% CALCULATE FRAME AVERAGE:

% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul = zeros(frame_Ysize,frame_Xsize);

for k = start_frame:end_frame  % loop through frames.
    % Get frame data:
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul = frame_accumul + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg = frame_accumul/(end_frame-start_frame+1);

imshow(frame_avg,[]); % display result image for graphical inspection.


%% SELECT REGION OF INTEREST:

% Select region of interest on frame average:
% user needs to place cursor markers on average image and input the
% following numbers into the program, eg. [1 512 256 512] to select the bottom half of image:
disp(' ') % empty line
limitsROI = input('Enter [xleft xright ytop ybottom] limits for region of interest: '); % request user input

xleft = limitsROI(1);
xright = limitsROI(2);
ytop = limitsROI(3);
ybottom = limitsROI(4);
% Save as parameters for output:
params.xleft_ROI = xleft;
params.xright_ROI = xright;
params.ytop_ROI = ytop;
params.ybottom_ROI = ybottom;

% frame average Region Of Interest:
frame_avgROI = frame_avg(ytop:ybottom,xleft:xright);

% SAVE current FIGURE (frame average) as a .png and then close the figure window:
figName = strcat('frame_Avg',image_label); % choose name of figure file for saving .png.
saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.

% Show region of interest only:
figure('position',[800 100 800 800]); % [left, bottom, width, height]; first two for lower-left corner of figure.
subplot(2,1,1)
imshow(frame_avgROI,[],'Border','tight');
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('frame_Avg_ROI',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.


%% SELECT BACKGROUND SPOT on ROI for background subtraction for all spots:

disp(' ') % empty line
bgnd_spot = input('Enter position [x_bgnd_spot y_bgnd_spot] for background spot within ROI: '); % request user input

x_bgnd_spot = bgnd_spot(1);
y_bgnd_spot = bgnd_spot(2);

% Save as parameters for output:
params.bgnd_spot_x = x_bgnd_spot;
params.bgnd_spot_y = y_bgnd_spot;



%% FIND CANDIDATE BRIGHT SPOTS FOR THE FRAME AVERAGE ROI:

% frame dimensions:
frame_Ysize = size(frame_avgROI,1);
frame_Xsize = size(frame_avgROI,2);
% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos (used in future sections):
[Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);
% Note that the image thresholding occurrs in two halves: separating top and bottom halves.
% Find candidate-bright-spots on first frame:
frame_Gray = mat2gray(frame_avgROI); % The input to function "findCandidateSpots" needs to be a grayscale image:

[candidate_spotsX_0 candidate_spotsY_0] = findCandidateSpots(frame_Gray,2); % Second input: use method 2, which seems to work better.
% See C:\Isabel\myMatlabFiles\findCandidateSpots.m.
% candidate_spotsX and candidate_spotsY are two column vectors of the same
% length containing the x and y coordinates of the candidate bright spots found on the image.
disp(' '); % empty line.
disp(['no. of candidate spots on frame average: ',num2str(length(candidate_spotsX_0))])

% Error control:
% Limit the max number of candidate spots (if eg. 260000 candidate spots are
% found, we will get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB").
% Select only the first params.max_num_candidates then.
if length(candidate_spotsX_0) > params.max_num_candidates
    candidate_spotsX_0 = candidate_spotsX_0(1:params.max_num_candidates);
    candidate_spotsY_0 = candidate_spotsY_0(1:params.max_num_candidates);
    disp(['NOTE!! no. of candidate spots has been limited to ',num2str(params.max_num_candidates)])
end
% Error control:
% If no candidates are found, break and exit:
if isempty(candidate_spotsX_0)
    disp('No candidates for bright spot founds. Exiting program.')
    return
end

% %Check graphically:
% figure;
% imshow(frame_Gray,[],'Border','tight');
% hold on;
% plot(candidate_spotsX_0,candidate_spotsY_0,'*');
% hold off;


%% FIND SPOT CENTRES AND DECIDE IF WE ACCEPT THEM OR NOT:

n =1; % Initialise index n (index for accepted spot centres which have a clipping flag equal to zero):
frame_to_search = frame_avgROI; % Initialise frame to search for spot centres.

disp(' '); % empty line.
disp('Finding spot centres... ');
disp(' '); % empty line.

% Find spot centre through iterative masking:
for m = 1:size(candidate_spotsX_0,1) % loop through all candidate spots.
    % Now find centre of bright spot using function findSpotCentre1frame:
    % use candidate spots as initial estimates and then iterate to find spot centre.
    % Image subarray ROI is a square of size 17x17 pixels (halfwidth is
    % 8 pixels), inner circular mask that moves inside the fixed 17x17
    % square has a radius of 5 pixels and the applied Gaussian
    % mask has a sigma of 2 pixels:
    spot_result = findSpotCentre1frame(frame_to_search,candidate_spotsX_0(m),candidate_spotsY_0(m),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,guess_sigma_Fit);
    % spot_result = findSpotCentre1frame(frame_to_search,candidate_spotsX_0(m),candidate_spotsY_0(m),2+params.square_halfsize,params.square_halfsize,params.square_halfsize/4);
    
    if (spot_result.ClipFlag == 0 && spot_result.noConverge == 0 && ...
            spot_result.SigmaFit <= params.sigmaFit_max && ...
            spot_result.SigmaFit >= params.sigmaFit_min && ...
            spot_result.SNR >= params.SNR_min &&...
            spot_result.rsqFit >= params.rsq_min)
        % Only accept and save result of found spot if clipping flag =0 and if values of sigmaFit, signal-to-noise ratio
        % and rsquare of fit are acceptable.
        spot_result.SpotNumber = n; % Add new field containing spot number to result structure.
        spot_final(n) = spot_result; % store "good" found spots in the final result spot_final, structure array.
        
        %-------------------------------
        if params.deflate==1 % see parameter section at the beginning.
            % "Deflation" process: subtract from raw frame image the corresponding
            % Gaussian fit of each found and accepted spot before finding next spot centre (enables acceptance of dimmer spots).
            %
            % Matrices containing the x and y positions in the image frame: Xpos and Ypos
            % Xpos is a matrix of the same size as frame, containing x values for all pixels and similarly for Ypos.
            % Calculate Xpos, Ypos at the beginning: [Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);
            % Parameters of Gaussian fit of previously accepted spot:
            x_fit = spot_final(n).CentreX;
            y_fit = spot_final(n).CentreY;
            I_fit = spot_final(n).I0Fit;
            sigma_fit = spot_final(n).SigmaFit;
            % deflated frame (frame with found spot subtracted):
            deflated_frame = frame_to_search - I_fit*exp(-((Xpos-x_fit).^2+(Ypos-y_fit).^2)/(2*sigma_fit^2));
            frame_to_search = deflated_frame; % update frame to search for finding next spot-centre.
            % % Graphical check of deflated frames:
            % subplot(1,2,1); imshow(frame_avgROI,[]); % frame is the original frame (always the same).
            % subplot(1,2,2); imshow(deflated_frame,[]);
        end
        %-------------------------------
        
        n = n+1; % advance index n for accepted spot centres.
    end
end

% % display the number of accepted spot-centres in frame average:
numAcceptedSpots = n-1;
disp(['no. of accepted spot centres in frame average: ',num2str(numAcceptedSpots)])
% Error control:
% If no spots are accepted, break and exit:
if numAcceptedSpots < 1
    disp(' ')
    disp('ERROR !!! ')
    disp('No spots accepted. Exiting program.')
    return
end


% convert results of found spot-centre positions to a useful form:
if (n-1) == 0 % error control: if no spots were accepted.
    found_spot_CentreX = [];
    found_spot_CentreY = [];
    % I need to create the whole spot_final structure with all its fields
    % here, just in case the number of accepted spots in the first frame is
    % zero, in order not to get error: "Subscripted assignment between
    % dissimilar structures".
    % Save empty spot (we need this, otherwise if in the last frame the no. of accepted spots is 0, there will be no result spot_final(end_frame,:) and the following functions will fail).
    spot_final(n).CentreX = [];
    spot_final(n).CentreY = [];
    spot_final(n).IspTot = [];
    spot_final(n).rsqFit = [];
    spot_final(n).SigmaFit = [];
    spot_final(n).I0Fit = [];
    spot_final(n).BgNoiseStd = [];
    spot_final(n).IbgAvg = [];
    spot_final(n).IbgTot = [];
    spot_final(n).SNR = [];
    spot_final(n).IinnerTot = [];
    spot_final(n).ClipFlag = [];
    spot_final(n).noConverge = [];
    spot_final(n).TrajNumber = [];
else
    found_spot_CentreX = [spot_final(:).CentreX]'; % column vector with found CentreX positions of all candidate spots.
    found_spot_CentreY = [spot_final(:).CentreY]'; % column vector with found CentreY positions of all candidate spots.
end

% Check graphically, plot accepted bright spots:
%figure;
subplot(2,1,2)
imshow(frame_Gray,[],'Border','tight');
hold on;
plot(found_spot_CentreX,found_spot_CentreY,'o','Color','g','MarkerSize',8) % plot accepted spot centres in green.
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% SAVE current FIGURE (frame avg and frame avg with spots overlaid) as a .png and then close the figure window:
figName = strcat('Avg_frame_Accepted_Spots',image_label); % choose name of figure file for saving .png.
saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('accepted_Bright_Spots',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% -----------------------------------------------------------------------


%% LOOK AT INDIVIDUAL FRAMES IN SEQUENCE
% Find integrated intensity from all spots.

% BACKGROUND CALCULATIONS:

% Loop through frames:
for k = start_frame:end_frame
    
    % Get frame data:
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    
    % Take ROI of each individual frame:
    frameROI = frame(ytop:ybottom,xleft:xright);
    
    % Find out BACKGROUND INTEGRATED INTENSITY:
    % Calculate sum of intensities in a square of size 2*params.square_halfsize+1
    % around selected bgnd spot position: x_bgnd_spot, y_bgnd_spot.
    % Create SQUARE mask of ones around bgnd spot.
    % "bgnd_mask" is of the same size as frameROI with
    % a small square of ones around the position of the bgnd spot and zeroes
    % elsewhere. Default half-size of small square is 4-5 pixels, given by "params.square_halfsize".
    bgnd_mask = zeros(size(frameROI,1),size(frameROI,2)); % pre-define image mask size, initialise.
    % Need error control when taking points at edge: to not use points too close to edge:
    if (round(x_bgnd_spot)-params.square_halfsize)<1 || (round(x_bgnd_spot)+params.square_halfsize)>size(frameROI,2) || (round(y_bgnd_spot)-params.square_halfsize)<1 || (round(y_bgnd_spot)+params.square_halfsize)>size(frameROI,1)
        % Do nothing (skip spot) if too close to edge. Matrix square_mask will remain zero.
    else
        for ii = (round(y_bgnd_spot)-params.square_halfsize):(round(y_bgnd_spot)+params.square_halfsize) % ii = index for rows
            for jj = (round(x_bgnd_spot)-params.square_halfsize):(round(x_bgnd_spot)+params.square_halfsize) % jj = index for columns
                bgnd_mask(ii,jj)=1;
            end
        end
        %     figure;
        %     imshow(bgnd_mask.*frameROI,[],'Border','tight')
        Ibgnd(k,1) = sum(sum(bgnd_mask.*frameROI)); % Ibgnd is a column vector with elements k=1:end_frame. Elements k<start_frame are zero. 
    end
end  
    
% Filter Ibgnd: 
[Ibgnd_filtered_0,tx,dx,sd,dsd,xpre] = chungKennedyFilter(Ibgnd(start_frame:end_frame),params.Wfilter,params.Rfilter);

% Normalise filtered-bgnd: correction factor to remove effects of varying
% excitation intensities with time:
for k = 1:length(Ibgnd_filtered_0)
    Ibgnd_filtered_normalised_0(k,1) = Ibgnd_filtered_0(k,1)/max(Ibgnd_filtered_0); % column vector, background normalised to its maximum value.
end

% Ibgnd_filtered_0 and Ibgnd_filtered_normalised_0 start at start_frame. They are column vectors of length (end_frame-start_frame+1).
% Pad with zeros at the beginning of vectors to make them of length end_frame.
Ibgnd_filtered = padarray(Ibgnd_filtered_0,[start_frame-1 0],'pre'); % prepend (start_frame-1) zeros to column vector.
Ibgnd_filtered_normalised = padarray(Ibgnd_filtered_normalised_0,[start_frame-1 0],'pre'); % prepend (start_frame-1) zeros to column vector.

% Save bgnd vector to params structure:
params.Ibgnd = Ibgnd; % column vector of length end_frame.
params.Ibgnd_filtered = Ibgnd_filtered; % starts at start_frame. Column vector of length (end_frame-start_frame+1).
params.Ibgnd_filtered_normalised = Ibgnd_filtered_normalised; % column vector of length (end_frame-start_frame+1).


% SPOT INTENSITY CALCULATIONS:

% Initialise vector of integrated spot intensities.
Ispots = zeros(end_frame,length(found_spot_CentreX)); 
% First index in Ispots (rows, moving vertically) is frame number (related to time).
% Second index is number of bright spot. So each column is the intensity
% frame by frame for a given spot number.
% Note that because we will then fill up this matrix starting at
% start_frame, the first few rows might be just zeroes.

% Loop through frames:
for k = start_frame:end_frame
    
    disp(['frame ',num2str(k)])
    
    % Get frame data:
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    
    % Take ROI of each individual frame:
    frameROI = frame(ytop:ybottom,xleft:xright);
    
    %     figure;
    %     imshow(frameROI,[],'Border','tight')
    
    % For each frame, loop through all accepteded spots:
    for n = 1:length(found_spot_CentreX)
        
        % disp(['spot number: ',num2str(n),' '])
        x_centre = found_spot_CentreX(n);
        y_centre = found_spot_CentreY(n);
        
        % Create SQUARE around position of found and accepted spots:
        % "square_mask" is a square image of the same size as frameROI with
        % a small square of ones around the position of the spot and zeroes
        % elsewhere. Default half-size of small square is 4-5 pixels, given by "params.square_halfsize".
        square_mask = zeros(size(frameROI,1),size(frameROI,2)); % pre-define image mask, square_mask, size, initialise.       
        % Need error control when taking points at edge: do not use points
        % too close to edge:
        if (round(x_centre)-params.square_halfsize)<1 || (round(x_centre)+params.square_halfsize)>size(frameROI,2) || (round(y_centre)-params.square_halfsize)<1 || (round(y_centre)+params.square_halfsize)>size(frameROI,1)
            % error control for when spots are too close to edge of image ROI.
            % Do nothing (skip spot) if too close to edge. Matrix square_mask will remain zero.
        else
            % only for a small square (saves time) of side 2*params.square_halfsize+1, around the
            % found spot, build square of ones:
            for ii = (round(y_centre)-params.square_halfsize):(round(y_centre)+params.square_halfsize) % ii = index for rows
                for jj = (round(x_centre)-params.square_halfsize):(round(x_centre)+params.square_halfsize) % jj = index for columns
                    square_mask(ii,jj)=1;
                end
            end
            %    nnz(square_mask) % number of non-zero elements in matrix square_mask.
            %     figure;
            %     imshow(square_mask,[],'Border','tight')
            %     figure;
            %     imshow(square_mask.*frameROI,[],'Border','tight')
            Ispots(k,n) = sum(sum(square_mask.*frameROI)); % integrated spot intensity.
            Ispots_bgndSubtracted(k,n) = Ispots(k,n)-Ibgnd_filtered(k); % subtract background (same for all spots).
            Ispots_corrected(k,n) = Ispots_bgndSubtracted(k,n)/Ibgnd_filtered_normalised(k); % correct for variations proportional to excitation intensity.
        end
    end
end

%------------------------------


%% Make result structure and eliminate empty spots (spots too close to figure edges which are all zeros):
% Final spots to analyse are saved in structure array "results".
% Also apply Chung-Kennedy filter to each trace and save into "results".

% 'params.Wfilter' and 'params.Rfilter' are the window size and weighting exponent of the
% Chung-Kennedy filter. See 'chungKennedyFilter.m'.
% Go to PARAMETERS section to find values of params.Wfilter and params.Rfilter.

nn = 1; % initialise final spot number (index nn).
for n=1:numAcceptedSpots % Loop through accepted spots (index n).
    if sum(Ispots(:,n))==0 % if all elements in column are equal to zero (add up to zero):
        % do nothing;
    else       
        results(nn).spot_number = nn;
        results(nn).Xcentre = found_spot_CentreX(n);
        results(nn).Ycentre = found_spot_CentreY(n);       
        results(nn).frame_number = (1:end_frame)'; % column vector
        results(nn).integrated_Ispot = Ispots(:,n); % column vector
        results(nn).integrated_Ispot_bgndSubtracted = Ispots_bgndSubtracted(:,n); % column vector
        results(nn).integrated_Ispot_corrected = Ispots_corrected(:,n); % column vector
        
        % I_to_use = results(nn).integrated_Ispot; % unfiltered data vector.
        I_to_use = results(nn).integrated_Ispot_corrected; % unfiltered data vector.
        % Use only data points starting at start_frame:
        [filtered_I_0,tx,dx,sd,dsd,xpre] = chungKennedyFilter(I_to_use(start_frame:end_frame),params.Wfilter,params.Rfilter);
        % Output 'filtered_I' is the filtered vector (ignore other outputs).
        % filtered_I starts at start_frame, it is a column vector of length (end_frame-start_frame+1).
        % Pad with zeros at the beginning to make it of length end_frame.
        filtered_I = padarray(filtered_I_0,[start_frame-1 0],'pre'); % prepend (start_frame-1) zeros to column vector.
        % Add fields to result structure:
        results(nn).filtered_Ispot = filtered_I; % starts at start_frame. Column vector of length end_frame.
        
        nn = nn + 1;
    end
end

%% Plot frame average ROI with numbers of final accepted spots (index nn) overlaid 
% on top to be able to visually identify them:
figure; 
imshow(frame_Gray,[],'Border','tight');
hold on;
for nn=1:length(results)
    text(results(nn).Xcentre,results(nn).Ycentre,num2str(nn),'Color',[1 1 0],'FontSize',8); % number in yellow.
end
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% SAVE current FIGURE as a .png and then close the figure window:
figName = strcat('Avg_frame_Accepted_Spot_Numbers',image_label); % choose name of figure file for saving .png.
saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.


%% Plot integrated intensity traces graphically as a 2D image, to get a quick view:

% results_Imatrix = [results.integrated_Ispot];
% figure;
% imshow(results_Imatrix,[min(results_Imatrix(results_Imatrix~=0)),max(results_Imatrix(results_Imatrix~=0))],'Border','tight','InitialMagnification',400);
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('IspotsMatrix_quickView',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% % -------------------------


%% Preliminary Graphical display: make figures with sets of 16 plots of integrated spot intensity
% versus frame number: 
nplots = ceil(length(results)/16);

for ifig=1:nplots
    figure('position',[500 100 1200 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
    for m=1:16
        l=16*(ifig-1)+ m;
        subplot(4,4,m); 
        if l<=length(results)
            plot(results(l).frame_number,results(l).integrated_Ispot_corrected,'-b');
            xlabel('frame no.');
            ylabel('Ispot tot');
            title(num2str(l));
            hold on;
            % plot(results(l).frame_number,results(l).integrated_Ispot_bgndSubtracted,'-r');
            % Plot Chung-Kennedy filtered (of corrected) trace:
            plot(results(l).frame_number,results(l).filtered_Ispot,'black-','LineWidth',1);
            % legend('bgnd subtracted','corrected','filtered','Location','Best');
            % legend('corrected','filtered','Location','Best');
            hold off;
            xlim([start_frame end_frame]);
        end
    end
    % Give time to user to look at figure:
    input('Look at plots for vertical limits for zooming later. Press any key to continue: '); 
    % SAVE current FIGURE as a .png and then close the figure window:
    figName = strcat('IspotVsTime_',image_label,'_',num2str(ifig)); % choose name of figure file for saving .png.
    saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
end
    
% Request user input for future Ylimits to save zoom of plots right after:
    limitsYzoom = input('Enter new vertical limits [ymin ymax] for zooming into plots: '); 

% Now make another figure zooming in/out, with those new Ylimits:
for ifig=1:nplots
    figure('position',[500 100 1200 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
    for m=1:16
        l=16*(ifig-1)+ m;
        subplot(4,4,m); 
        if l<=length(results)
            plot(results(l).frame_number,results(l).integrated_Ispot_corrected,'-b');
            xlabel('frame no.');
            ylabel('Ispot tot');
            title(num2str(l));
            hold on; 
            % plot(results(l).frame_number,results(l).integrated_Ispot_bgndSubtracted,'-r');
            % Plot Chung-Kennedy filtered (of corrected) trace:
            plot(results(l).frame_number,results(l).filtered_Ispot,'black-','LineWidth',1);
            % legend('bgnd subtracted','corrected','filtered','Location','Best');
            % legend('corrected','filtered','Location','Best');
            hold off;
            xlim([start_frame end_frame]);
            ylim([limitsYzoom(1) limitsYzoom(2)]); % zoom in vertically using limitsYzoom entered by user.
        end
    end
    % SAVE current FIGURE as a .png and then close the figure window:
    figName = strcat('IspotVsTime_zoom_',image_label,'_',num2str(ifig)); % choose name of figure file for saving .png.
    saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
end


%% LOOP THROUGH ALL SPOTS: loop through all Ispot-versus-frame traces:
% Calculate distribution of pair-wise differences and save into "results".
% Calculate Fourier transform of the former and detect peaks in it.

% Loop through final accepted spots (index nn) and apply filter:
for nn=1:length(results)

    I_to_use = results(nn).integrated_Ispot_corrected; % unfiltered data vector of length end_frame.
    % Begin large figure for a given spot:
    % figure('position',[200 100 1500 1000]); 
    figure('position',[10 10 2200 1100]); % [left, bottom, width, height]; first two for lower-left corner of figure.
   
    % subplot 1:
    subplot(3,5,1) % Integrated Ispot vs frame number with Chung Kennedy filter, bgnd subtracted and corrected::
    plot(results(nn).frame_number,I_to_use,results(nn).frame_number,results(nn).filtered_Ispot,'-k','LineWidth',1);
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),' bgnd subtracted & corrected']});
    xlim([start_frame end_frame]);
    % subplot 2:
    subplot(3,5,2) % Zoom in/out, integrated Ispot vs frame number with Chung Kennedy filter:
    plot(results(nn).frame_number,I_to_use,results(nn).frame_number,results(nn).filtered_Ispot,'-k','LineWidth',1);
    xlim([start_frame end_frame]);
    ylim([limitsYzoom(1) limitsYzoom(2)]); % zoom in vertically using limitsYzoom entered by user.
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),' bgnd subtracted & corrected']});      
    
    % subplot 3:
    subplot(3,5,3) % Zoom in/out, uncorrected Ispot vs frame number.
    plot(results(nn).frame_number,results(nn).integrated_Ispot_bgndSubtracted,'-r');
    xlim([start_frame end_frame]);
    ylim([limitsYzoom(1) limitsYzoom(2)]); % zoom in vertically using limitsYzoom entered by user.
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),' bgnd subtracted, not corrected']});  
    
    % subplot 4:
    subplot(3,5,4) % Plot integrated background intensity versus frame number:
    plot(results(nn).frame_number,params.Ibgnd,results(nn).frame_number,params.Ibgnd_filtered,'-k'); 
    v = axis; % v contains [xmin xmax ymin ymax] of zoomed-out plot.
    xlim([start_frame end_frame]);
    ylim([v(3) v(4)]); % use limits of zoom out, otherwise, changing the xlimits in previous line automatically zooms it in "y" too.
    xlabel('frame no.'); ylabel('Integrated Background Intensity');    
            
    % Calculate and plot distribution of pair-wise differences (see fullPwD.m):
    % points_per_bin = 20; See PARAMETERS section.
    % Use only data between start_frame and end_frame.
    [binCentres freqCounts] = fullPwD(I_to_use(start_frame:end_frame),params.nbins,0); % Calculate distribution of pair-wise differences.
    % Add fields to result structure:
    results(nn).binCentres_Ipwd = binCentres'; % column vector.
    results(nn).freqCounts_Ipwd = freqCounts'; % column vector.
    % Plot distributions of pair-wise differences:
    % subplot 5:
    subplot(3,5,5) % full histogram.
    bar(binCentres,freqCounts,'r'); % plot a bar graph of the full histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title({['spot no. ',num2str(nn)]});
    ylim([0 1.1*max(freqCounts)]); % re-scale vertical axis.
    
    pos_negatives = find(binCentres < 0); % positions of negative bin-centre values.
    pos_positives = find(binCentres > 0); % positions of positive bin-centre values.
    if isempty(pos_negatives) % this is a fast dodgy way of sorting this out for the moment.
        pos_negatives = 1;
        warning('No negative side of histogram!');
    end
    if isempty(pos_positives) % this is a fast dodgy way of sorting this out for the moment.
        pos_positives = 1;
        warning('No positive side of histogram!');
    end
    binCentres_neg = binCentres(pos_negatives);
    freqCounts_neg = freqCounts(pos_negatives);
    binCentres_pos = binCentres(pos_positives);
    freqCounts_pos = freqCounts(pos_positives);
    
    % subplot 6:
    subplot(3,5,6) % :
    bar(binCentres_neg,freqCounts_neg,'b'); % plot negative histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title({['spot no. ',num2str(nn),' prob distrib of I drops']});
    
    % subplot 7:
    subplot(3,5,7) % :
    bar(binCentres_pos,freqCounts_pos,'g'); % plot positive histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title({['spot no. ',num2str(nn),' prob distrib of I jumps > 0']});
     
    % Calculate and plot power spectrum (ps) of distribution of pair-wise
    % differences of each spot, and peaks found in it, in order of
    % increasing spectral power (see FourierAndFindPeaks.m):
    % params.x_limit_spectr = -0.7*min(binCentres); % this is a good estimate of the max value we want to cut the spectrum at in the x axis.
    [ps_x ps_y ps_peaks_x ps_peaks_y] = FourierAndFindPeaks(binCentres,freqCounts,0,0);
    % last input = 0 means figure is not displayed within function "FourierAndFindPeaks.m".
    % Add fields to result structure:
    results(nn).power_spectrum_X = ps_x; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y = ps_y; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X = ps_peaks_x; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y = ps_peaks_y; % column vector: y-values of peaks found in power spectrum (spectral power).

    % Plot spectral power versus intensity_step values:
    % subplot 8:
    subplot(3,5,8) % Power spectrum of prob distrib. of I-pair-wise differences:
    plot(ps_x,ps_y,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of FULL prob distrib. of I-pair-wise differences');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps_peaks_x),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps_peaks_x(i),3));
        text(ps_peaks_x(i),ps_peaks_y(i),label_to_print,'FontSize',10);
    end
   
    % Now select and use ONLY NEGATIVE values (intensity drops) of pairwise intensity differences (horiz
    % axis in the previous histogram):
    % Calculate and plot power spectrum (ps) of negative side of distribution of pair-wise
    % differences for each spot, and peaks found in it (see FourierAndFindPeaks.m):
    [ps_x_neg ps_y_neg ps_peaks_x_neg ps_peaks_y_neg] = FourierAndFindPeaks(binCentres_neg,freqCounts_neg,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum_X_neg = ps_x_neg; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y_neg = ps_y_neg; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X_neg = ps_peaks_x_neg; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y_neg = ps_peaks_y_neg; % column vector: y-values of peaks found in power spectrum (spectral power).
    
    % Plot spectral power versus intensity_step values for negative side only:
    % subplot 9:
    subplot(3,5,9) % Power spectrum of prob distrib. of I-pair-wise differences < 0:
    plot(ps_x_neg,ps_y_neg,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of prob distrib. of I-pair-wise-differences < 0');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow_neg = min(length(ps_peaks_x_neg),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow_neg
        label_to_print_neg = strcat('\leftarrow',num2str(ps_peaks_x_neg(i),3));
        text(ps_peaks_x_neg(i),ps_peaks_y_neg(i),label_to_print_neg,'FontSize',10);
    end   
%     axis_limits_current_figure = axis; % gives limits of axis of current figure: [xmin xmax ymin ymax].
%     text(0.1*axis_limits_current_figure(2),0.9*axis_limits_current_figure
%     (4),{['spot number ',num2str(nn)]}); % write at top left of figure the spot number.


    % Now select and use ONLY POSITIVE values (intensity jumps>0) of pairwise intensity differences:
    % Calculate and plot power spectrum (ps) of positive side of distribution of pair-wise
    % differences for each spot, and peaks found in it:
    [ps_x_pos ps_y_pos ps_peaks_x_pos ps_peaks_y_pos] = FourierAndFindPeaks(binCentres_pos,freqCounts_pos,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum_X_pos = ps_x_pos; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y_pos = ps_y_pos; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X_pos = ps_peaks_x_pos; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y_pos = ps_peaks_y_pos; % column vector: y-values of peaks found in power spectrum (spectral power).
    
    % Plot spectral power versus intensity_step values for positive side only:
    % subplot 10:
    subplot(3,5,10) % Power spectrum of prob distrib. of I-pair-wise differences > 0:
    plot(ps_x_pos,ps_y_pos,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of prob distrib. of I-pair-wise-differences > 0');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow_pos = min(length(ps_peaks_x_pos),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow_pos
        label_to_print_pos = strcat('\leftarrow',num2str(ps_peaks_x_pos(i),3));
        text(ps_peaks_x_pos(i),ps_peaks_y_pos(i),label_to_print_pos,'FontSize',10);
    end   

    % METHOD 2: pair-wise differences of pair-wise differences.
    PD = allPD(I_to_use(start_frame:end_frame)); % Get column vector with all Intensity Pair-wise Differences (PD). See allPD.m.
    PD_neg = PD(PD < 0); % negative pair-wise differences.
    PD_pos = PD(PD > 0); % positive pair-wise differences.
    if isempty(PD_neg) 
        warning('No negative pair-wise differences found!');
    end
    if isempty(PD_pos) 
        warning('No positive pair-wise differences found!');
    end

    % Histogram of pair-wise differences of NEGATIVE pair-wise differences:
    [binCentres2_neg freqCounts2_neg] = fullPwD(PD_neg,params.nbins,0); % Calculate distribution of pair-wise differences of the negative pair-wise differences.
    % Add fields to result structure:
    results(nn).binCentres2_neg = binCentres2_neg'; % column vector.
    results(nn).freqCounts2_neg = freqCounts2_neg'; % column vector.
    % subplot 11:
    subplot(3,5,11) % Histogram of pair-wise differences of negative pair-wise differences:
    bar(binCentres2_neg,freqCounts2_neg,'b'); % plot a bar graph of the histogram.
    xlabel('pair-wise differences of intensity pair-wise differences < 0 ');
    ylabel('frequency');
    title({['spot no. ',num2str(nn)]});
    ylim([0 1.1*max(freqCounts2_neg)]); % re-scale vertical axis.
       
    % Calculate and plot power spectrum (ps2 for method 2) of previous histogram:
    [ps2_x_neg ps2_y_neg ps2_peaks_x_neg ps2_peaks_y_neg] = FourierAndFindPeaks(binCentres2_neg,freqCounts2_neg,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum2_X_neg = ps2_x_neg; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum2_Y_neg = ps2_y_neg; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum2_peaks_X_neg = ps2_peaks_x_neg; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum2_peaks_Y_neg = ps2_peaks_y_neg; % column vector: y-values of peaks found in power spectrum (spectral power).
    % Plot spectral power versus intensity_step values:
    % subplot 12:
    subplot(3,5,12) % Power spectrum of prob distrib. of pair-wise differences of negative Intensity pair-wise differences:
    plot(ps2_x_neg,ps2_y_neg,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of pair-wise diffs of I-pair-wise diffs <0');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps2_peaks_x_neg),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps2_peaks_x_neg(i),3));
        text(ps2_peaks_x_neg(i),ps2_peaks_y_neg(i),label_to_print,'FontSize',10);
    end
    
    % Histogram of pair-wise differences of POSITIVE pair-wise differences:
    [binCentres2_pos freqCounts2_pos] = fullPwD(PD_pos,params.nbins,0); % Calculate distribution of pair-wise differences of the positive pair-wise differences.
    % Add fields to result structure:
    results(nn).binCentres2_pos = binCentres2_pos'; % column vector.
    results(nn).freqCounts2_pos = freqCounts2_pos'; % column vector.
    % subplot 13:
    subplot(3,5,13) % Histogram of pair-wise differences of positive pair-wise differences:
    bar(binCentres2_pos,freqCounts2_pos,'b'); % plot a bar graph of the histogram.
    xlabel('pair-wise differences of intensity pair-wise differences > 0');
    ylabel('frequency');
    title({['spot no. ',num2str(nn)]});
    ylim([0 1.1*max(freqCounts2_pos)]); % re-scale vertical axis.   
    
    % Calculate and plot power spectrum (ps2 for method 2) of previous histogram:
    [ps2_x_pos ps2_y_pos ps2_peaks_x_pos ps2_peaks_y_pos] = FourierAndFindPeaks(binCentres2_pos,freqCounts2_pos,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum2_X_pos = ps2_x_pos; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum2_Y_pos = ps2_y_pos; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum2_peaks_X_pos = ps2_peaks_x_pos; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum2_peaks_Y_pos = ps2_peaks_y_pos; % column vector: y-values of peaks found in power spectrum (spectral power).
    % Plot spectral power versus intensity_step values:
    % subplot 14:
    subplot(3,5,14) % Power spectrum of prob distrib. of pair-wise differences of positive intensity pair-wise differences:
    plot(ps2_x_pos,ps2_y_pos,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of pair-wise diffs of I-pair-wise diffs >0');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps2_peaks_x_pos),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps2_peaks_x_pos(i),3));
        text(ps2_peaks_x_pos(i),ps2_peaks_y_pos(i),label_to_print,'FontSize',10);
    end
    
    
    % Histogram of pair-wise differences of ALL pair-wise differences:
    [binCentres2 freqCounts2] = fullPwD(PD,params.nbins,0); % Calculate distribution of pair-wise differences of all pair-wise differences.
    % Add fields to result structure:
    results(nn).binCentres2 = binCentres2'; % column vector.
    results(nn).freqCounts2 = freqCounts2'; % column vector.  
    
    % Calculate and plot power spectrum (ps2 for method 2) of previous histogram:
    [ps2_x ps2_y ps2_peaks_x ps2_peaks_y] = FourierAndFindPeaks(binCentres2,freqCounts2,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum2_X = ps2_x; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum2_Y = ps2_y; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum2_peaks_X = ps2_peaks_x; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum2_peaks_Y = ps2_peaks_y; % column vector: y-values of peaks found in power spectrum (spectral power).
    % Plot spectral power versus intensity_step values:
    % subplot 15:
    subplot(3,5,15) % Power spectrum of prob distrib. of pair-wise differences of positive Intensity pair-wise differences:
    plot(ps2_x,ps2_y,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    title('Power Spectrum of pair-wise diffs of ALL I-pair-wise diffs');
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps2_peaks_x),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps2_peaks_x(i),3));
        text(ps2_peaks_x(i),ps2_peaks_y(i),label_to_print,'FontSize',10);
    end
    
    % SAVE current FIGURE as a .png and then close the figure window:
    figName = strcat('results_',image_label,'_spot_',num2str(nn)); % choose name of figure file for saving .png.
    saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.

end
% -------------------------


%% Get SUM of ALL integrated INTENSITIES from all spots for each frame and analyse it:
% % plot vs time to try and get photobleaching steps: sum over all accepted Ispots...
% 
% % This gives a column vector in which each element corresponds to a frame and is the sum of all
% % integrated spot intensities (from final accepted spots) in that frame:
% sumIspot_all_spots = sum(results_Imatrix,2); % vector, sum across second dimension of matrix.
% 
% % Apply Chung-Kennedy filter to the sum trace:
% % Use only data points starting at start_frame.
% [filtered_Isum,tx,dx,sd,dsd,xpre] = chungKennedyFilter(sumIspot_all_spots(start_frame:end_frame),params.Wfilter,params.Rfilter);
% % Output 'filtered_Isum' is the filtered vector (ignore other outputs).
% frame_num_for_Isum_filtered = (start_frame:end_frame)';
% 
% 
% figure;
% plot(sumIspot_all_spots);
% xlabel('frame no.');
% ylabel('sum Ispot all spots');
% title('sum of integrated intensities of all accepted spots');
% % Plot Chung-Kennedy filtered trace:
% hold on;
% plot(frame_num_for_Isum_filtered,filtered_Isum,'black-','LineWidth',2)
% hold off;
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('I_allSpots',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% 
% 
% % Full calculation of Pair-wise Differences (PwD) for column vector
% % "sumIspot_all_spots", calculation of their probability distribution 
% % and optional graph display (last input param = 1) of their histograms:
% [binCentres_all_spots freqCounts_all_spots] = fullPwD(sumIspot_all_spots,params.nbins,1);
% % See fullPwD.m :
% % [binCentres freqCounts] = fullPwD(input_vector,params.nbins,display_figures)
% 
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('Ipwd_all_spots',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% 
% % Note that here all pair-wise differences are calculated for all pairs, not only for
% % consecutive pairs of points.
% 
% % binCentres: row vector containing the resulting bin centres of the full
% % histogram of all pair-wise differences.
% % freqCounts: row vector containing the occurrences (or frequency) for each
% % value of pair-wise differences.
% % The first plot shows the full histogram.
% % The second plot shows the highest end of the full histogram, related to a
% % possible large background level of the initial signal in input_vector,
% % when no background subtraction has been applied.
% % The third plot shows the lower end of the full histogram, related to
% % photobleaching steps, blinking, noise, etc.
% 
% % Calculate and plot power spectrum (ps) of distribution of pair-wise
% % differences of all spots together, and peaks found in it:
% x_limit = -0.5*min(binCentres_all_spots); % this is a good estimate of the value we want to cut the spectrum at in the x axis.
% [ps_x_all_spots ps_y_all_spots ps_peaks_x_all_spots ps_peaks_y_all_spots] = FourierAndFindPeaks(binCentres_all_spots,freqCounts_all_spots,1,x_limit);
% % last input = 1 means figure is displayed.
% 
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('power_spectrum_allSpots',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% 
% 
% % -------------------------
% % Now select and use only negative values (intensity drops) of pairwise intensity differences (horiz
% % axis in the previous histogram):
% pos_negatives = find(binCentres_all_spots < 0); % positions of negative bin-centre values.
% if isempty(pos_negatives) % this is a fast dodgy way of sorting this out for the moment.
%     pos_negatives = 1;
%     warning('No negative side of histogram!');
% end
% binCentres_all_spots_neg = binCentres_all_spots(pos_negatives);
% freqCounts_all_spots_neg = freqCounts_all_spots(pos_negatives);
% % Calculate and plot power spectrum (ps) of distribution of pair-wise
% % differences of all spots together for negative side of spectrum only, and peaks found in it:
% x_limit = -0.5*min(binCentres_all_spots_neg); % this is a good estimate of the value we want to cut the spectrum at in the x axis.
% [ps_x_all_spots_neg ps_y_all_spots_neg ps_peaks_x_all_spots_neg ps_peaks_y_all_spots_neg] = FourierAndFindPeaks(binCentres_all_spots_neg,freqCounts_all_spots_neg,1,x_limit);
% % last input = 1 means figure is displayed.
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('power_spectrum_allSpots_neg',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% 
% % ---------------------

%% INTERPOLATE power spectra to have them all equally uniformly sampled: 
% % This in order to be able to take their average:
% 
% xAll = [results.power_spectrum_X]; % matrix where each column is the x-axis for the power spectrum for each spot.
% x1min = min(min(xAll)); % Choose the overall minimum of all columns.
% x1max = max(xAll(2,:)); % The maximum in all columns is Infinity, so choose next value...
% x2 = (x1min:1:x1max)'; % new x-axis, column vector, uniformly sampled x-axis, it is the same for all spots.
% % idem for spectrum of only negative side of distribution of pair-wise
% % differences. (Need to do it slightly different because x-axes can be of different lengths).
% xAll_neg = []; % initialise vector column before loop.
% for nn = 1:length(results)
%     xAll_neg = [xAll_neg ; results(nn).power_spectrum_X_neg]; % make a column vector which adds one by one the column vectors power_spectrum_X_neg for all spots.
% end
% x1min_neg = min(xAll_neg); % Choose the overall minimum of all.
% x1max_neg = max(xAll_neg(xAll_neg~=Inf)); % The maximum in all columns is Infinity, so remove Inf before finding max.
% x2_neg = (x1min_neg:1:x1max_neg)'; % new x-axis, column vector, uniformly sampled x-axis, it is the same for all spots.
% 
% for nn = 1:length(results)
%     x1_0 = results(nn).power_spectrum_X; % column vector, original x-axis (non-uniformly sampled) of power spectrum.
%     y1_0 = results(nn).power_spectrum_Y; % column vector, original y-axis of power spectrum.
%     x1 = x1_0(2:end); % eliminate first point, which is Inf.
%     y1 = y1_0(2:end);
%     % Interpolate: the new interpolated spectrum is x2, y2:
%     y2 = interp1(x1,y1,x2); % interpolated y-axis: one point for each value of x2.
%     % Note that y-values outside the range of x1 are returned as Nan in y2,
%     % but we can later operate with vectors ok and take their mean despite this.
%     % Add fields to result structure:
%     results(nn).interpol_spectr_X = x2; % column vector: x-axis of interpolated power spectrum (Istep values).
%     results(nn).interpol_spectr_Y = y2; % column vector: y-axis of interpolated power spectrum (spectral power).
%     
%     % idem for power spectrum of "negative" part:
%     x1_0_neg = results(nn).power_spectrum_X_neg; % column vector, original x-axis (non-uniformly sampled) of power spectrum.
%     y1_0_neg = results(nn).power_spectrum_Y_neg; % column vector, original y-axis of power spectrum.
%     x1_neg = x1_0_neg(2:end); % eliminate first point, which is Inf.
%     y1_neg = y1_0_neg(2:end);
%     % Interpolate: the new interpolated spectrum is x2, y2:
%     y2_neg = interp1(x1_neg,y1_neg,x2_neg); % interpolated y-axis: one point for each value of x2.
%     % Add fields to result structure:
%     results(nn).interpol_spectr_X_neg = x2_neg; % column vector: x-axis of interpolated power spectrum for negative part (Istep values).
%     results(nn).interpol_spectr_Y_neg = y2_neg; % column vector: y-axis of interpolated power spectrum for negative part (spectral power).
% end


%% Calculate and display AVERAGE spectrum by taking the mean of all spectra from all spots: 

% % The mean spectrum averaging for all spots is avg_spectr_X, avg_spectr_Y:
% spectr_Matrix = [results.interpol_spectr_Y]; % Matrix where each column vector is the y-data (spectral power) for a given spot.
% avg_spectr_Y = mean(spectr_Matrix,2); % column vector, averaged y-axis of the mean spectrum.
% avg_spectr_X = results(1).interpol_spectr_X; % column vector, averaged x-axis of the mean spectrum (the x-axis is the same for all spots).
% 
% % Idem for the mean spectrum, averaging for all spots, of negative part: avg_spectr_X_neg, avg_spectr_Y_neg:
% spectr_Matrix_neg = [results.interpol_spectr_Y_neg]; % Matrix where each column vector is the y-data (spectral power) for a given spot.
% avg_spectr_Y_neg = mean(spectr_Matrix_neg,2); % column vector, averaged y-axis of the mean spectrum.
% avg_spectr_X_neg = results(1).interpol_spectr_X_neg; % column vector, averaged x-axis of the mean spectrum (the x-axis is the same for all spots).
% 
% % Display and save graphs:
% figure('position',[1200 100 500 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.;
% subplot(2,1,1)
% plot(avg_spectr_X,avg_spectr_Y); 
% xlim([0 1000]);
% subplot(2,1,2)
% plot(avg_spectr_X_neg,avg_spectr_Y_neg);
% xlim([0 1000]);
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('Avg_power_spectrum',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.


%% Put all Intensity pair-wise differences (PD) from all spots together in one vector,
% % calculate the probability distribution of that and its Fourier transform:
% % Label "together" for all figures in this section.
% 
% % initialise vector column before loop:
% PD_all_spots_together = []; % all intensity pair-wise differences (PDs) together.
% for nn = 1:length(results)
%     PD_spot = allPD(results(nn).integrated_Ispot); % column vector with PDs for spot number "nn".
%     PD_all_spots_together = [PD_all_spots_together; PD_spot]; % make a column vector which adds one by one the column vectors power_spectrum_X_neg for all spots.
% end
% 
% % Calculate probability distribution of all intensity pair-wise differences
% % (PD) from all spots together and plot as a histogram:
% [freqCounts_together binCentres_together] = hist(PD_all_spots_together,300); % Get histogram. Second input is number of bins.
%  
% figure('position',[1200 100 500 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
% subplot(3,1,1) % full histogram.
% bar(binCentres_together,freqCounts_together,'r'); % plot a bar graph of the full histogram.
% xlabel('Intensity pair-wise differences all spots TOGETHER'); ylabel('frequency');
% ylim([0 max(freqCounts_together)]); % re-scale vertical axis.
% 
% % Separate components in previous histogram by changing range to plot:
% subplot(3,1,2) % Background component (higher end of histogram):
% bar(binCentres_together,freqCounts_together,'b'); % plot same full histogram.
% xlabel('Intensity pair-wise differences TOGETHER');
% ylabel('frequency');
% title('background component');
% xlim([0.3*max(binCentres_together) 1.1*max(binCentres_together)]); % re-scale horizontal axis.
% 
% subplot(3,1,3) % Photobleaching/blinking/noise component (lower end of histogram):
% bar(binCentres_together,freqCounts_together,'g'); % plot same full histogram.
% xlabel('Intensity pair-wise differences TOGETHER');
% ylabel('frequency');
% title('photobleaching/blinking/noise component');
% xlim([min(binCentres_together) -min(binCentres_together)]); % re-scale horizontal axis.
% 
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('Ipwd_histogram_together',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
% 
% % Replot same histogram and zoom in vertically and horizontally:
% figure;
% bar(binCentres_together,freqCounts_together,'g'); % plot same full histogram.
% xlabel('Intensity pair-wise differences');
% ylabel('frequency');
% title('photobleaching/blinking/noise component');
% xlim([min(binCentres_together) 0]); % re-scale horizontal axis.
% ylim([0 200]);
% 
% 
% % Calculate and plot power spectrum (ps) of distribution of pair-wise
% % differences and peaks found in it:
% x_limit = -0.5*min(binCentres_together); % this is a good estimate of the value we want to cut the spectrum at in the x axis.
% [ps_x_together ps_y_together ps_peaks_x_together ps_peaks_y_together] = FourierAndFindPeaks(binCentres_together,freqCounts_together,1,x_limit);
% % last input = 1 means figure is displayed.
% 
% % SAVE current FIGURE as a .png and then close the figure window:
% figName = strcat('power_spectrum_together',image_label); % choose name of figure file for saving .png.
% saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.


%% FINAL OUTPUT RESULT:

% Final result:

result_final = {results params};


%% ------------------------------------
% Determine on-off bliking times.

