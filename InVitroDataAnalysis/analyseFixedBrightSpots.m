function result_final = analyseFixedBrightSpots(image_label,start_frame,end_frame,folder_label)
%
% March 2012. Created by Isabel Llorente Garcia. 
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function has a lot in common with findTrajects.m
%
% This function is used to analyse fluorescent labels which are fixed
% (glued) onto a glass surface when imaged. This means tracking is not necessary.
% I first get a frame average and find the position of the fixed bright
% spots and then monitor the fluorescence in those spots over time.
% The aim is to extract the typical intensity of a single-molecule of a
% given fluorophore.
% The background subtraction is local, around the position of the bright
% spots found in the frame average, using a circular region and a slightly bigger square
% region around it.
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
%     power_spectrum_X_pos
%     power_spectrum_Y_pos
%     power_spectrum_peaks_X_pos
%     power_spectrum_peaks_Y_pos
%
% "params": is a structure containing as fields the different parameters
% used for the analysis:
%                     image_label: 
%                    folder_label: 
%                     start_frame: 
%                       end_frame: 
%              max_num_candidates: 
%                    sigmaFit_min: 
%                    sigmaFit_max: 
%                         SNR_min: 
%                         rsq_min: 
%                         deflate: 
%                         Wfilter: 
%                         Rfilter: 
%                           nbins: 
%                  x_limit_spectr: 
%          subtract_bgnd_fit_hist: 
%                    use_filtered:
%                      image_path: string
%                       numFrames: 
%                     frame_Ysize: 
%                     frame_Xsize: 
%                       xleft_ROI: 
%                      xright_ROI: 
%                        ytop_ROI: 
%                     ybottom_ROI: 
%                     bgnd_spot_x: 
%                     bgnd_spot_y: 
%                           Ibgnd: vector
%                binCentres_Ibgnd: vector
%                freqCounts_Ibgnd: vector
%           power_spectrum_X_bgnd: vector
%           power_spectrum_Y_bgnd: vector
%     power_spectrum_peaks_X_bgnd: vector
%     power_spectrum_peaks_Y_bgnd: vector
%                  Ibgnd_rmsNoise: 
%                 Ibgnd_ampli_fit: 
    
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
% These parameters are for the Gaussian masking method to find the
% centre of the spot, and they also affect the total integrated spot
% intensity calculation and background subtraction. 

% PARAMETERS for deciding if we accept a spot centre found by the function
% findSpotCentre1frame or not:
params.sigmaFit_min = -inner_circle_radius; % minimum acceptable sigma of gaussian fit to spot, in pixels (2).
params.sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4).
params.SNR_min = 2; % minimum acceptable signal-to-noise ratio (at least 2) (SNR as defined in findSpotCentre1frame.m).
params.rsq_min = 0.2; % minimum acceptable r-square value (0.2) (goodness of gaussian fit to spot).

% PARAMETER deflate: to decide if deflation method is applied or not (subtract each
% found spot to allow detection of weaker spots):
params.deflate = 1; % 1 for yes, 0 for no.

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
params.x_limit_spectr = 2000; % (5000, 10000 or 3000 for Rodrigo's data, 1200-3000 for Oxphos data GFP, 600 mCherry)

% PARAMETER: subtract (1) or not (0) Gaussian fit to distrib of Ibackground
% pair-wise differences (both to plot histograms and to calculate their spectrum and find peaks):
params.subtract_bgnd_fit_hist = 0;

% PARAMETER: use filtered intensity trace (and background-corrected) for
% analysis:
params.use_filtered = 0;

% PARAMETER: do intensity pair wise differences (1) for analysis or not (0):
params.doPairWiseDiffs = 0;

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
frame_accumul = zeros(frame_Ysize,frame_Xsize); % this is of class double.

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
% % Request user input:
limitsROI = input('Enter [xleft xright ytop ybottom] limits for region of interest: '); % request user input
% % or) Used a fixed ROI, no user input requested:
% limitsROI = [250 500 250 500];

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


% %% SELECT BACKGROUND SPOT on ROI for background subtraction for all spots:
% This is not used anymore... (Isabel, December 2012).
% disp(' ') % empty line
% bgnd_spot = input('Enter position [x_bgnd_spot y_bgnd_spot] for background spot within ROI: '); % request user input
% 
% x_bgnd_spot = bgnd_spot(1);
% y_bgnd_spot = bgnd_spot(2);
% 
% % Save as parameters for output:
% params.bgnd_spot_x = x_bgnd_spot;
% params.bgnd_spot_y = y_bgnd_spot;



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
    % Image subarray ROI is a square of size (2*subarray_halfwidth)^2 pixels (halfwidth is
    % subarray_halfwidth pixels), inner circular mask that moves inside the
    % fixed square has a radius of 5 pixels and the applied Gaussian
    % mask has a sigma of inner_circle_radius pixels:
    spot_result = findSpotCentre1frame(frame_to_search,candidate_spotsX_0(m),candidate_spotsY_0(m),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,guess_sigma_Fit);
    
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

% ---------------------------
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
        % elsewhere. Default half-size of small square is 4-5 pixels, given by "params.subarray_halfwidth".
        square_mask = zeros(size(frameROI,1),size(frameROI,2)); % pre-define image mask, square_mask, size, initialise.       
        % Create smaller CIRCULAR mask of radius "params.inner_circle_radius":
        % "circle_mask" is a square image of the same size as frameROI with
        % a small circle of ones around the position of the spot and zeroes
        % elsewhere:
        circle_mask = zeros(size(frameROI,1),size(frameROI,2)); % pre-define circular mask, initialise;
        % Need error control when taking points at edge: do not use points
        % too close to edge:
        if (round(x_centre)-params.subarray_halfwidth)<1 || (round(x_centre)+params.subarray_halfwidth)>size(frameROI,2) || (round(y_centre)-params.subarray_halfwidth)<1 || (round(y_centre)+params.subarray_halfwidth)>size(frameROI,1)
            % error control for when spots are too close to edge of image ROI.
            % Do nothing (skip spot) if too close to edge. Matrix square_mask will remain zero.
        else
            % only for a small square (saves time) of side 2*params.subarray_halfwidth+1, around the
            % found spot, build square of ones:
            for ii = (round(y_centre)-params.subarray_halfwidth):(round(y_centre)+params.subarray_halfwidth) % ii = index for rows
                for jj = (round(x_centre)-params.subarray_halfwidth):(round(x_centre)+params.subarray_halfwidth) % jj = index for columns
                    % squared subarray mask around spot centre:
                    square_mask(ii,jj)=1;
                    % inner-circle mask: if distance to point centre is <=inner_circle_radius, then 1, otherwise, 0:
                    sum_sq = (jj-x_centre)^2 + (ii-y_centre)^2;
                    if round(sqrt(sum_sq))<= params.inner_circle_radius
                        circle_mask(ii,jj)=1;
                    else
                        circle_mask(ii,jj)=0;
                    end
                end
            end
            num_pixels_inCircle(k,n) = nnz(circle_mask); % number of pixels within circular spot-signal mask (number of non-zero elements in matrix).
            
            % Define background mask: subtract circle_mask from square_mask
            % to be left with ones only at the edges of the square, outside
            % the inner circle:
            bgnd_mask = square_mask - circle_mask;      
            num_pixels_inBgnd(k,n) = nnz(bgnd_mask); % number of pixels within background mask (number of non-zero elements in matrix).
            %             figure;
            %             imshow(bgnd_mask.*frameROI,[],'Border','tight')
            pos_bgnd = find(bgnd_mask==1); % positions of bgnd only intensities.
            I1 = bgnd_mask.*frameROI; % bgnd region, image.
            Ibg_tot(k,n) = sum(I1(pos_bgnd)); % Sum of all intenstities in bgnd region. Ibgnd is a column vector with elements k=1:end_frame. Elements k<start_frame are zero.
            Ibg_avg(k,n) = mean(I1(pos_bgnd)); % Average background intensity per pixel, use mean for the bgnd to exclude hot pixels or nearby bright spots.
            Ibgnd(k,n) = Ibg_avg(k,n)*num_pixels_inCircle(k,n); % Equivalent total background intensity in spot region.
            
            Ispots(k,n) = sum(sum(circle_mask.*frameROI)); % total spot intensity before bgnd subtraction (Iinner_total).
            
            % Calculate background-corrected ROI image:
            I2 = frameROI-Ibg_avg(k,n);
            % Calculate standard deviation of remaining noise in background
            % corrected image I2:
            bg_noise_offset_afterBGsubtract(k,n) = mean(I2(pos_bgnd)); % offset noise level per pixel after background subtraction, should be close to zero.
            bg_noise_std(k,n) = std(I2(pos_bgnd)); % standard deviation of matrix elements in bgnd region.
            % Total spot intensity (within the inner circle mask), background corrected:
            Ispots_bgndSubtracted(k,n) = sum(sum(I2.*circle_mask)); % integrated spot intensity.
            
            %    nnz(square_mask) % number of non-zero elements in matrix square_mask.
            %     figure;
            %     imshow(square_mask,[],'Border','tight')
            %     figure;
            %     imshow(square_mask.*frameROI,[],'Border','tight')   
            
%             figure; imshow(square_mask,[]); figure; imshow(circle_mask,[]); figure; imshow(bgnd_mask,[]);
%             figure; imshow(bgnd_mask.*frameROI,[]); figure; imshow(circle_mask.*frameROI,[]); figure; imshow(I2.*circle_mask,[]);
        end
    end
end

% -------------------------------------



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
        results(nn).frame_number = (start_frame:end_frame)'; % column vector        
        results(nn).integrated_Ispot = Ispots(start_frame:end_frame,n); % column vector, total spot intensity before bgnd subtraction.
        results(nn).integrated_Ispot_bgndSubtracted = Ispots_bgndSubtracted(start_frame:end_frame,n); % column vector, total spot intensity after bgnd subtraction.
        results(nn).Ibgnd = Ibgnd(start_frame:end_frame,n); % column vector, equivalent total bgnd intensity within spot region/circular mask.
        results(nn).Ibg_tot = Ibg_tot(start_frame:end_frame,n); % column vector, total intensity in bgnd mask region.
        results(nn).Ibg_avg = Ibg_avg(start_frame:end_frame,n); % column vector, average intensity per pixel in bgnd mask region.
        results(nn).num_pixels_inCircle = num_pixels_inCircle(start_frame:end_frame,n); % column vector, number of pixels in inner circle region, spot signal region.
        results(nn).num_pixels_inBgnd = num_pixels_inBgnd(start_frame:end_frame,n); % column vector, number of pixels in bgnd mask around spot.
        results(nn).bg_noise_offset_afterBGsubtract = bg_noise_offset_afterBGsubtract(start_frame:end_frame,n); % column vector, this is ~0.
        results(nn).bg_noise_std = bg_noise_std(start_frame:end_frame,n); % column vector, bgnd noise standard deviation within bgnd region.
        % note: all vectors of length end_frame-start_frame+1.
        
        % ------------
        % BACKGROUND CALCULATIONS:      
        Ibgnd_vector = results(nn).Ibgnd; % column vector, equivalent total bgnd intensity within spot region/circular mask.
          
        if params.doPairWiseDiffs == 0
            % Calculate histogram of intensity trace (minus its average value) directly, without doing
            % pair-wise differences, so histogram of intensity levels:
            [freqCounts_bgnd,binCentres_bgnd] = hist(Ibgnd_vector-mean(Ibgnd_vector),params.nbins);            
        else
            % Calculate distribution of pair-wise differences of background intensity, Ibgnd (see fullPwD.m):
            % Use only data between start_frame and end_frame.
            [binCentres_bgnd freqCounts_bgnd] = fullPwD(Ibgnd_vector,params.nbins,0); % Calculate distribution of pair-wise differences.       
        end
      
        % Calculate power spectrum (ps) of distribution of pair-wise differences (or of intensity levels) of Ibgnd:
        [ps_x_bgnd ps_y_bgnd ps_peaks_x_bgnd ps_peaks_y_bgnd] = FourierAndFindPeaks(binCentres_bgnd,freqCounts_bgnd,0,0);
        % last input = 0 means figure is not displayed within function "FourierAndFindPeaks.m".
        
        % Fit distrib of Ibgnd pair-wise differences/intensity levels to a Gaussian (Poissonian shot noise with large photon no.):
        fun_for_fit = fittype('ampli*exp(-I^2/(2*sigma^2))','independent','I'); % define Gaussian funtion to fit to, with 'I' as independent variable;
        options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
        options.StartPoint = [max(freqCounts_bgnd) std(Ibgnd_vector)]; % give guess parameters for fit. This avoids a warning message. Guess for I0 is max(data_to_fit) and guess for sigma is standard deviation of original data.
        % options.Lower = [min_sigma_fit min_ampli_fit]; % Lower bounds for fit parameters (see PARAMETER section at beginning of this file).
        % options.Upper = [max_sigma_fit max_ampli_fit]; % Upper bounds for fit parameters (see PARAMETER section at beginning of this file).
        [fit_result gof] = fit(binCentres_bgnd',freqCounts_bgnd',fun_for_fit,options); % do fit. fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the good of fitness.
        % fit_param_names = coeffnames(fit_result); % get fit parameter names to know their order: first one is I0, second one is sigma_spot.
        fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is "ampli", second one is "sigma".
        ampli_bgnd_fit = fit_param_values(1); % fitted amplitude from fit.
        sigma_bgnd_fit = fit_param_values(2); % Sigma of Gaussian fit, i.e. approx rms size (stdev) of background intensity noise.
        % rsq_fit = gof.rsquare; % rsquare coefficient of fit.
        % errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
        % errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
        
        % Save bgnd vector and other results to params structure within results:
        results(nn).binCentres_Ibgnd = binCentres_bgnd'; % column vector.
        results(nn).freqCounts_Ibgnd = freqCounts_bgnd'; % column vector.
        results(nn).power_spectrum_X_bgnd = ps_x_bgnd; % column vector: x-axis of power spectrum (Istep values).
        results(nn).power_spectrum_Y_bgnd = ps_y_bgnd; % column vector: y-axis of power spectrum (spectral power).
        results(nn).power_spectrum_peaks_X_bgnd = ps_peaks_x_bgnd; % column vector: x-values of peaks found in power spectrum (Istep values).
        results(nn).power_spectrum_peaks_Y_bgnd = ps_peaks_y_bgnd; % column vector: y-values of peaks found in power spectrum (spectral power).
        results(nn).Ibgnd_rmsNoise = sigma_bgnd_fit;
        results(nn).Ibgnd_ampli_fit = ampli_bgnd_fit;
        % ------------
        
        
        % I_to_use = results(nn).integrated_Ispot; % data vector.
        I_to_use = results(nn).integrated_Ispot_bgndSubtracted; % data vector.
        % Use only data points starting at start_frame:
        [filtered_I_0,tx,dx,sd,dsd,xpre] = chungKennedyFilter(I_to_use,params.Wfilter,params.Rfilter);
        % Output 'filtered_I' is the filtered vector (ignore other outputs).
        % filtered_I starts at start_frame, it is a column vector of length (end_frame-start_frame+1).
        filtered_I = filtered_I_0;
%         % Pad with zeros at the beginning to make it of length end_frame.
%         filtered_I = padarray(filtered_I_0,[start_frame-1 0],'pre'); % prepend (start_frame-1) zeros to column vector.
        % Add fields to result structure:
        results(nn).filtered_Ispot = filtered_I; % starts at start_frame. Column vector of length end_frame-start-frame+1.
        
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
            plot(results(l).frame_number,results(l).integrated_Ispot_bgndSubtracted,'-b');
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
    
%     input('Look at plots for vertical limits for zooming later. Press any key to continue: '); 
    
    % SAVE current FIGURE as a .png and then close the figure window:
    figName = strcat('IspotVsTime_',image_label,'_',num2str(ifig)); % choose name of figure file for saving .png.
    saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
end
    
% % Request user input for future Ylimits to save zoom of plots right after:
limitsYzoom = input('Enter new vertical limits [ymin ymax] for zooming into plots: '); 
    
% % Use fixed limits, no user input requested:
% limitsYzoom = [-2000 6000];
    

% Now make another figure zooming in/out, with those new Ylimits:
for ifig=1:nplots
    figure('position',[500 100 1200 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
    for m=1:16
        l=16*(ifig-1)+ m;
        subplot(4,4,m); 
        if l<=length(results)
            plot(results(l).frame_number,results(l).integrated_Ispot_bgndSubtracted,'-b');
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
% Calculate distribution of pair-wise differences or of intensity levels and save into "results".
% Calculate Fourier transform of the former and detect peaks in it.

% Loop through final accepted spots (index nn) and apply filter:
for nn=1:length(results)

    if params.use_filtered == 1
        I_to_use = results(nn).filtered_Ispot; % filtered data vector.
    else
        I_to_use = results(nn).integrated_Ispot_bgndSubtracted; % unfiltered data vector of length end_frame.
    end
    
    % Begin large figure for a given spot:
    % figure('position',[200 100 1500 1000]); 
    figure('position',[10 10 2200 1100]); % [left, bottom, width, height]; first two for lower-left corner of figure.
   
    % subplot 1:
    subplot(3,4,1) % Integrated Ispot vs frame number with Chung Kennedy filter, bgnd subtracted and corrected::
    plot(results(nn).frame_number,I_to_use,results(nn).frame_number,results(nn).filtered_Ispot,'-k','LineWidth',1);
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),' bgnd subtracted']});
    xlim([start_frame end_frame]);
    % subplot 2:
    subplot(3,4,2) % Zoom in/out, integrated Ispot vs frame number with Chung Kennedy filter:
    plot(results(nn).frame_number,I_to_use,results(nn).frame_number,results(nn).filtered_Ispot,'-k','LineWidth',1);
    xlim([start_frame end_frame]);
    ylim([limitsYzoom(1) limitsYzoom(2)]); % zoom in vertically using limitsYzoom entered by user.
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),' bgnd subtracted']});      
    
    % subplot 3:
    subplot(3,4,3) % Zoom in/out, uncorrected Ispot vs frame number.
    plot(results(nn).frame_number,results(nn).integrated_Ispot,'-r');
    xlim([start_frame end_frame]);
    xlabel('frame no.'); ylabel('Integrated spot Intensity');
    title({['spot no. ',num2str(nn),'before bgnd subtraction']});  
    
    % subplot 4:
    subplot(3,4,4) % Plot integrated background intensity versus frame number:
    plot(results(nn).frame_number,results(nn).Ibgnd,'-k'); 
    v = axis; % v contains [xmin xmax ymin ymax] of zoomed-out plot.
    xlim([start_frame end_frame]);
    ylim([0 v(4)]); % use limits of zoom out, otherwise, changing the xlimits in previous line automatically zooms it in "y" too.
    xlabel('frame no.'); ylabel('Integrated mean background intensity');    
    title('I background equivalent');
            
    if params.doPairWiseDiffs == 0
        % Calculate histogram of intensity trace directly, without doing
        % pair-wise differences, so histogram of intensity levels:
        [freqCounts,binCentres] = hist(I_to_use,params.nbins);
        xlabel_string = 'Intensities';
        % Add fields to result structure:
        results(nn).binCentres_Ilevels = binCentres'; % column vector.
        results(nn).freqCounts_Ilevels = freqCounts'; % column vector.
    else
        % Calculate and plot distribution of pair-wise differences (see fullPwD.m):
        % points_per_bin = 20; See PARAMETERS section.
        % Use only data between start_frame and end_frame.
        [binCentres freqCounts] = fullPwD(I_to_use,params.nbins,0); % Calculate distribution of pair-wise differences.
        xlabel_string = 'Intensity pair-wise differences';
        % Add fields to result structure:
        results(nn).binCentres_Ipwd = binCentres'; % column vector.
        results(nn).freqCounts_Ipwd = freqCounts'; % column vector.
    end
    
    % Separate positive and negative pair-wise differences/intensity levels: 
    pos_negatives = find(binCentres < 0); % positions of negative bin-centre values.
    pos_positives = find(binCentres > 0); % positions of positive bin-centre values.
    if isempty(pos_negatives) % this is a fast dodgy way of sorting this out for the moment.
        % pos_negatives = 1;
        disp('!!!  No negative side of histogram !!!');
        binCentres_neg = [0];
        freqCounts_neg = [0];
    else
        binCentres_neg = binCentres(pos_negatives);
        freqCounts_neg = freqCounts(pos_negatives);
    end
    
    if isempty(pos_positives) % this is a fast dodgy way of sorting this out for the moment.
        % pos_positives = 1;
        disp('!!! No positive side of histogram !!!');
        binCentres_pos = [0];
        freqCounts_pos = [0];
    else
        binCentres_pos = binCentres(pos_positives);
        freqCounts_pos = freqCounts(pos_positives);
    end
    
    % Plot all pair-wise differences/intensity levels:
    if params.subtract_bgnd_fit_hist == 1
        % Subtract Gaussian fit to distribution of background intensities:
        y_bgnd_GaussianFit = ampli_bgnd_fit*exp(-binCentres.^2/(2*sigma_bgnd_fit^2)); % Gaussian fitted shape of histogram for background, to subtract.
        freqCounts_bgndFitSubtracted = freqCounts - y_bgnd_GaussianFit;
        freqCounts_to_use = freqCounts_bgndFitSubtracted;
    else
        freqCounts_to_use = freqCounts;
    end
    % Plot distributions of pair-wise differences/intensity levels:
    % subplot 5:
    subplot(3,4,5) % full histogram.
    bar(binCentres,freqCounts_to_use,'r'); % plot a bar graph of the full histogram.
    xlabel(xlabel_string);
    ylabel('frequency');
    title({['spot no. ',num2str(nn)]});
    
    
    % Plot only negative pair-wise differences/intensity levels:
    if params.subtract_bgnd_fit_hist == 1
        % Subtract Gaussian fit to distribution of background intensities:
        y_bgnd_GaussianFit_neg = ampli_bgnd_fit*exp(-binCentres_neg.^2/(2*sigma_bgnd_fit^2)); % Gaussian fitted shape of histogram for background, to subtract.
        freqCounts_bgndFitSubtracted_neg = freqCounts_neg - y_bgnd_GaussianFit_neg;
        freqCounts_neg_to_use = freqCounts_bgndFitSubtracted_neg; % reasign value
    else
        freqCounts_neg_to_use = freqCounts_neg;
    end
    % subplot 6:
    subplot(3,4,6) % :
    bar(binCentres_neg,freqCounts_neg_to_use,'b'); % plot negative histogram.
    xlabel(xlabel_string);
    ylabel('frequency');
    if params.doPairWiseDiffs == 0
        title({['spot no. ',num2str(nn),' prob distrib of I < 0']});
    else
        title({['spot no. ',num2str(nn),' prob distrib of I drops']});
    end
        
    % Plot only positive pair-wise differences/intensity levels:
    if params.subtract_bgnd_fit_hist == 1
        % Subtract Gaussian fit to distribution of background intensities:
        y_bgnd_GaussianFit_pos = ampli_bgnd_fit*exp(-binCentres_pos.^2/(2*sigma_bgnd_fit^2)); % Gaussian fitted shape of histogram for background, to subtract.
        freqCounts_bgndFitSubtracted_pos = freqCounts_pos - y_bgnd_GaussianFit_pos;
        freqCounts_pos_to_use = freqCounts_bgndFitSubtracted_pos; % reasign value
    else
        freqCounts_pos_to_use = freqCounts_pos;
    end
    % subplot 7:
    subplot(3,4,7) % :
    bar(binCentres_pos,freqCounts_pos_to_use,'g'); % plot positive histogram.
    xlabel(xlabel_string);
    ylabel('frequency');
    if params.doPairWiseDiffs == 0
        title({['spot no. ',num2str(nn),' prob distrib of I > 0']});
    else
        title({['spot no. ',num2str(nn),' prob distrib of I jumps > 0']});
    end
     
    % Calculate and plot power spectrum (ps) of distribution of pair-wise
    % differences/intensity levels of each spot, and peaks found in it, in order of
    % increasing spectral power (see FourierAndFindPeaks.m):
    % params.x_limit_spectr = -0.7*min(binCentres); % this is a good estimate of the max value we want to cut the spectrum at in the x axis.
    [ps_x ps_y ps_peaks_x ps_peaks_y] = FourierAndFindPeaks(binCentres,freqCounts_to_use,0,0);
    % last input = 0 means figure is not displayed within function "FourierAndFindPeaks.m".
    % Add fields to result structure:
    results(nn).power_spectrum_X = ps_x; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y = ps_y; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X = ps_peaks_x; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y = ps_peaks_y; % column vector: y-values of peaks found in power spectrum (spectral power).
   
    % Plot spectral power versus intensity_step values:
    % subplot 9:
    subplot(3,4,9) % Power spectrum of prob distribution of ALL I-pair-wise differences/intensity levels:
    plot(ps_x,ps_y,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    if params.doPairWiseDiffs == 0
        title('Power Spectrum of prob distrib. of intensities');
    else
        title('Power Spectrum of prob distrib. of ALL I-pair-wise diffs');
    end
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps_peaks_x),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps_peaks_x(i),3));
        text(ps_peaks_x(i),ps_peaks_y(i),label_to_print,'FontSize',10);
    end
    
    % Now select and use ONLY NEGATIVE values (intensity drops) of pairwise intensity differences/intensity levels (horiz
    % axis in the previous histogram):
    % Calculate and plot power spectrum (ps) of negative side of distribution of pair-wise
    % differences for each spot, and peaks found in it (see FourierAndFindPeaks.m):
    [ps_x_neg ps_y_neg ps_peaks_x_neg ps_peaks_y_neg] = FourierAndFindPeaks(binCentres_neg,freqCounts_neg_to_use,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum_X_neg = ps_x_neg; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y_neg = ps_y_neg; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X_neg = ps_peaks_x_neg; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y_neg = ps_peaks_y_neg; % column vector: y-values of peaks found in power spectrum (spectral power).
    
    % Plot spectral power versus intensity_step values for negative side only:
    % subplot 10:
    subplot(3,4,10) % Power spectrum of prob distrib. of I-pair-wise differences < 0:
    plot(ps_x_neg,ps_y_neg,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    if params.doPairWiseDiffs == 0
        title('Power Spectrum of prob distrib. of I < 0');
    else
        title('Power Spectrum of prob distrib. of I-pair-wise-differences < 0');
    end
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


    % Now select and use ONLY POSITIVE values (intensity jumps>0) of pairwise intensity differences/intensity levels:
    % Calculate and plot power spectrum (ps) of positive side of distribution of pair-wise
    % differences for each spot, and peaks found in it:
    [ps_x_pos ps_y_pos ps_peaks_x_pos ps_peaks_y_pos] = FourierAndFindPeaks(binCentres_pos,freqCounts_pos_to_use,0,0);
    % Add fields to result structure:
    results(nn).power_spectrum_X_pos = ps_x_pos; % column vector: x-axis of power spectrum (Istep values).
    results(nn).power_spectrum_Y_pos = ps_y_pos; % column vector: y-axis of power spectrum (spectral power).
    results(nn).power_spectrum_peaks_X_pos = ps_peaks_x_pos; % column vector: x-values of peaks found in power spectrum (Istep values).
    results(nn).power_spectrum_peaks_Y_pos = ps_peaks_y_pos; % column vector: y-values of peaks found in power spectrum (spectral power).
    
    % Plot spectral power versus intensity_step values for positive side only:
    % subplot 11:
    subplot(3,4,11) % Power spectrum of prob distrib. of I-pair-wise differences > 0:
    plot(ps_x_pos,ps_y_pos,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    if params.doPairWiseDiffs == 0
        title('Power Spectrum of prob distrib. of I > 0');
    else
        title('Power Spectrum of prob distrib. of I-pair-wise-differences > 0');
    end
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow_pos = min(length(ps_peaks_x_pos),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow_pos
        label_to_print_pos = strcat('\leftarrow',num2str(ps_peaks_x_pos(i),3));
        text(ps_peaks_x_pos(i),ps_peaks_y_pos(i),label_to_print_pos,'FontSize',10);
    end          
   
    % Plot distributions of pair-wise differences or of intensity levels of Ibgnd:
    % subplot 8:
    subplot(3,4,8) % full histogram.
    bar(results(nn).binCentres_Ibgnd,results(nn).freqCounts_Ibgnd,'c'); % plot a bar graph of the full histogram.
    hold on;
    y_axis_GaussianFit = results(nn).Ibgnd_ampli_fit*exp(-results(nn).binCentres_Ibgnd.^2/(2*(results(nn).Ibgnd_rmsNoise)^2));
    plot(results(nn).binCentres_Ibgnd,y_axis_GaussianFit,'-k'); % Superimpose Gaussian fit.
    hold off;
    if params.doPairWiseDiffs == 0
        xlabel('Intensity level');
        title({['spot no. ',num2str(nn),' Ibgnd-mean(Ibgnd) ']});
    else
        xlabel('Ibgnd pair-wise differences');
        title({['spot no. ',num2str(nn),' Ibgnd ']});
    end
    ylabel('frequency');
    ylim([0 1.1*max(results(nn).freqCounts_Ibgnd)]); % re-scale vertical axis.
    xlim([-4*results(nn).Ibgnd_rmsNoise 4*results(nn).Ibgnd_rmsNoise]);  % re-scale horizontal axis.    
    
    % Plot power spectrum (ps) of distribution of pair-wise differences of Ibgnd:
    % Plot spectral power versus Ibgnd_step values:
    % subplot 12:
    subplot(3,4,12) % Power spectrum of prob distrib. of ALL Ibgnd pair-wise differences:
    plot(ps_x_bgnd,ps_y_bgnd,'b-','LineWidth',1.5);
    xlim([0 params.x_limit_spectr]);
    if params.doPairWiseDiffs == 0
        title('Power Spectrum of prob distrib. of all Ibgnd intensities');
    else
        title('Power spectrum of prob distrib. of ALL Ibgnd-pair-wise diffs');
    end
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps_peaks_x_bgnd),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps_peaks_x_bgnd(i),3));
        text(ps_peaks_x_bgnd(i),ps_peaks_y_bgnd(i),label_to_print,'FontSize',10);
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

