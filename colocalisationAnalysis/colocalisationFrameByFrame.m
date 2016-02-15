function result_struct = colocalisationFrameByFrame(image_label,bf_image_label,start_frame,end_frame,start_frame_forFrameAvg,end_frame_forFrameAvg,roi,output_label,colocalis_threshold,factor_forThresholding)
%
% Isabel Llorente-Garcia: August 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Analysis frame by frame of possible colocalisation or
% anti-colocalisation of bright spots in the top and bottom channels of
% fluorescence images.
% See
% C:\Isabel\DataAnalysis\myAnalysis\OxphosColocalisationAnalysis.doc
%
% Inputs:
% - image_label (fluorescence image sequence) such as '513', '490', etc... which corresponds to a certain .sif number of image sequence
% - bf_image_label: bright field image sequence (usually 10 frames only).
% - start_frame,end_frame: frames to use for the frame average. 
% - start_frame_forFrameAvg and end_frame_forFrameAvg: frames used to
% calculate the frame average which is used to generate the cell-signal and
% background masks.
% "end_frame_forFrameAvg" can be set to 'end' when all frames in the image sequence, from start_frame to
% 'end', want to be used.
% - roi: region of interest: region aroun the cell on
% either the top or bottom channels (referred to an image which is of size
% of only the top or bottom half). "roi" is a row vector [xleft xright ytop ybottom].
% Be slightly generous with the roi around the cell, to accommodate for the
% shift between top and bottom channels (could be up to ~10 pixels in
% either direction).
% Enter 'full_image' as "roi" if you want to use the full image size as
% region of interest for both channels.
% - output_label: string to label output .png files.
% - colocalis_threshold (0.144): values between -+ colocalis_threshold are
% considered colocalised. Positive final "colocalisation values" mean bright stuff in top channel, negative ones
% mean bright stuff in bottom channel, and values around zero mean bright
% stuff in both channels, ie., colocalised.
% - factor_forThresholding (1.5): parameter for thresholding the image to get the cell-signal mask within
% function getCellMaskAndBoundary3.m. The larger it its, the higher the
% threshold and the smaller the cell-signal mask region.
% If equal to 1, function getCellMaskAndBoundary3.m is equivalent to
% function getCellMaskAndBoundary.m
% 
%
% Example of how to call this function:
% rs498 = colocalisationFrameByFrame('498','497',5,30,5,500,'full_image',' ') for fluorescence image sequence 498, finding the displacement between top and bottom channels from bright field image '497'. .
% Another example: rs498_A = colocalisationFrameByFrame('498','497',5,30,5,'end',[10 60 10 60],'A'). Using
% last frame as end_frame_forFrameAvg.
% rs498 = colocalisationFrameByFrame('498','497',7,7+25,7,'end',[307 384 37 129],'')
%
% Output: 
% Graphic .png files are saved with results from the colocalisation
% analysis for each frame (colocalisation graph, histogram of
% colocalisation values and colocalisation fractions).
% it is a result structure with as many elements as frames analysed
% (start_frame to end_frame). Each element has fields:
% image_label
% bf_image_label
% start_frame
% end_frame
% start_frame_forFrameAvg
% end_frame_forFrameAvg
% roi
% output_label
% translateBy
% frameNumber
% colocalis_values % vector.
% colocalis_threshold
% colocalised_fraction
% topChannel_fraction
% bottomChannel_fraction


%% Create new directory for saving plots and results:

% Create new folder for outputs inside current folder:
output_folder_name = 'colocalisResults';
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(output_folder_name); % make new directory.


%% Calculate shift between top and bottom channels using bright field image-sequence: 
 
% Use the result of function overlapBF.m here. Calculates the shift needed
% for best overlap of the top and bottom channels of the bright field (bf)
% image, using all frames (from 1 to 'end') of the bright field image
% sequence.
% Reminder:  [shift_y shift_x] =
% overlapBF(image_label,start_frame,end_frame,exclude_y_top,exclude_y_bot,exclude_x_left,exclude_x_right,show_plots).
disp('Bright field image:')
[y1 x1] = overlapBF(bf_image_label,1,'end',30,30,30,30,0,0); % excluding 30 pixels on each edge. (30 isabel)
translateBy = [y1 x1];
% ---------------------
% % Or choose by hand the amount by which to shift the image (down and left):
% translateBy = [8,-1];
% % ---------------------
% % If you want to skip the translation:
% translateBy = [0 0];
% % ---------------------
disp('Top channel shifted by :');
disp(translateBy)




%% Operations on frame average (average of frames between start_frame_forFrameAvg and end_frame_forFrameAvg):
% ----------------------------------------------------------------------

% 1 - Calculate frame average between start_frame_forFrameAvg and end_frame_forFrameAvg:

disp('Fluorescence image:')
fr_avg = frameAverage(image_label,start_frame_forFrameAvg,end_frame_forFrameAvg,0,0); 


% 2 -  Separate top and bottom channels (usually red and green channels):

% Middle row in image that separates top half (red channel) from bottom
% half (green channel). If the image is 512x512, half = 256.
half = round(size(fr_avg,1)/2); % round to nearest integer: takes care of sizes in pixels which are odd numbers.
whole = size(fr_avg,1); % full vertical size.
% Separate top and bottom channels (usually red and green channels)
% and convert to grayscale (0 to 1) images (so final result is only
% qualitative):
Btop0 = fr_avg(1:half,:); % top channel.
Bbot0 = fr_avg((half+1):whole,:); % bottom channel.
% "0" at end means full image before taking the roi.


% 3 -  Shift the top image:

% Shift the top channel to counteract slight displacement between channels:
% eg. [5,-4]: shift down 5 pixels and left 4 pixels.
% Shifted top image:
Btop02 = circshift(Btop0,translateBy); % shift top image by vector translateBy.
% "0" at end means full image before taking the roi.
result1 = Btop02-Bbot0; % shifted top image minus bottom one.
% result1 has values between -1 and 1.
% Bright areas with values close to 1 will correspond to spots found in the
% top frame but not in the bottom one. 
% Dark areas with values close to -1 will correspond to spots found in the
% bottom frame but not in the top one.
% Grey areas with values close to 0 will correspond to either the
% background or areas where bright spots are colocalised in both channels.


% 4 - Use a region of interest on the frame average:
% The region of interest is given by the input "roi".

xleft = roi(1);
xright = roi(2);
ytop = roi(3);
ybottom = roi(4);

if strcmp(roi,'full_image') % if input "roi" is 'full_image'
    Btop2 = Btop02;
    Bbot = Bbot0;
else % if input "roi" is 'full_image' a row vector [xleft xright ytop ybottom]
    Btop2 = Btop02(ytop:ybottom,xleft:xright);
    Bbot = Bbot0(ytop:ybottom,xleft:xright);
end
% I keep the "2" at the end of Btop to indicate that it is the shifted top channel.


% 5 - Re-scale frame-average images to make them between 0 and 1:

% Note that the re-scaling is done after shifting and taking the roi.
Btop2 = mat2gray(Btop2); % convert to grayscale.
Bbot = mat2gray(Bbot); % convert to grayscale.


% 6 - Get signal and background masks (cell foreground shapes):

[SignalMaskTop2 cellBoundaryTop2]= getCellMaskAndBoundary3(Btop2,factor_forThresholding); % for shifted top channel, after taking roi.
[SignalMaskBot cellBoundaryBot]= getCellMaskAndBoundary3(Bbot,factor_forThresholding); % for bottom channel, after taking roi.
BgndMaskTop2 = ~SignalMaskTop2; %  background mask is the complementary of the signal mask, shifted top channel, after taking roi.
BgndMaskBot = ~SignalMaskBot; % background mask is the complementary of the signal mask, bottom channel, after taking roi.

full_signalMask = (SignalMaskTop2 | SignalMaskBot); % union of signal masks of top and bottom channels.
full_bgndMask = (BgndMaskTop2 & BgndMaskBot); % intersection of background masks of top and bottom channels.

% % Alternative, to test threshold with single-label strains: use signal mask
% % of only top channel and not its union with the bottom mask:
% full_signalMask = SignalMaskTop2;
% full_bgndMask = BgndMaskTop2;


%% Read in the image-sequence data:
%----------------------------------

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.



%% Operations frame by frame (loop between start_frame and end_frame):
%--------------------------------------------------------------------

i = 1; % index for saving to result structure, "result_structs".
all_colocal_values = []; % initialise empty vector to accumulate colocalisation values.

for n = start_frame:end_frame % loop through chosen frames
    
    % 1 - Get single frame data:
    frame = image_data(n).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    disp(strcat('Frame number: ',num2str(n))) 
    
    % 2 -  Separate top and bottom channels (usually red and green
    % channels). Use previous values of "half" and "whole".
    % Separate top and bottom channels (usually red and green channels)
    % and convert to grayscale (0 to 1) images (so final result is only
    % qualitative):
    Btop0_singleFrame = frame(1:half,:); % top channel.
    Bbot0_singleFrame = frame((half+1):whole,:); % bottom channel.
    % "0" at end means full image before taking the roi.
    
    % 3 -  Shift the top image:
    % Shift the top channel to counteract slight displacement between channels:
    % eg. [5,-4]: shift down 5 pixels and left 4 pixels.
    % Shifted top image:
    Btop02_singleFrame = circshift(Btop0_singleFrame,translateBy); % shift top image by vector translateBy.
    % "0" at end means full image before taking the roi.
    result1_singleFrame = Btop02_singleFrame-Bbot0_singleFrame; % shifted top image minus bottom one.
    
    % 4 - Use a region of interest on the frame average:
    % The region of interest is given by the input "roi".
    % Use previous values of xleft,xright,ytop,ybottom.
    if strcmp(roi,'full_image') % if input "roi" is 'full_image'
        Btop2_singleFrame = Btop02_singleFrame;
        Bbot_singleFrame = Bbot0_singleFrame;
    else % if input "roi" is 'full_image' a row vector [xleft xright ytop ybottom]
        Btop2_singleFrame = Btop02_singleFrame(ytop:ybottom,xleft:xright);
        Bbot_singleFrame = Bbot0_singleFrame(ytop:ybottom,xleft:xright);
    end
    % I keep the "2" at the end of Btop to indicate that it is the shifted top channel.
    
    % 5 - Re-scale frame-average images to make them between 0 and 1:
    % Note that the re-scaling is done after shifting and taking the roi.
    Btop2_singleFrame = mat2gray(Btop2_singleFrame); % convert to grayscale.
    Bbot_singleFrame = mat2gray(Bbot_singleFrame); % convert to grayscale. 
       
    
    % Colocalisation result: 
    % ----------------------
    % Use the cell-signal and background masks obtained from the frame average, from the re-scaled roi
    % cuts of top and bottom channels in the frame average;
    result2_prelim = Btop2_singleFrame-Bbot_singleFrame; % preliminary result.
    max1 = max(max(result2_prelim));
    min1 = min(min(result2_prelim));
    
    % "result2" below is the final result, for a colour plot. The background is now
    % lower than any values in the foreground (-2), values close to +1 (or
    % actually max1) correspond to spots found in the top frame but not in the bottom one,
    % values close to -1 (or actually min1) correspond to spots found in the
    % bottom frame but not the top one, and values close to 0 correspond to
    % co-localised spots.
    % Note: (BgndMaskTop2 & BgndMaskBot) is the final background mask: intersection of the translated top one and the bottom one (ones in
    % bgnd and zeros in cell/foreground region).
    % (SignalMaskTop2 | SignalMaskBot) is the union of the two signal cell
    % masks.
    result2 = full_signalMask.*(Btop2_singleFrame-Bbot_singleFrame) - 2.*full_bgndMask;
    result2_topPlot = full_signalMask.*Btop2_singleFrame - full_bgndMask;
    result2_botPlot = full_signalMask.*Bbot_singleFrame - full_bgndMask;
    
    figure;
    subplot(3,2,1); imshow(result2_topPlot,[0,1]); colorbar('location','eastoutside'); title('shifted top channel');
    subplot(3,2,3); imshow(result2_botPlot,[0,1]); colorbar('location','eastoutside'); title('bottom channel');
    subplot(3,2,5); imshow(result2,[-max(abs(min1),max1),max(abs(min1),max1)]);
    colorbar('location','eastoutside'); title('shifted top minus bottom');
    colormap(jet);
     
    
    % Output, colocalisation result:
    colocalis_values = result2(full_signalMask==1); % column vector. Select only values inside full signal mask. Colocalisation values (between -1 and 1).
    all_colocal_values = [all_colocal_values; colocalis_values]; % accumulate values, append new column vector.
    
    % Plot histogram: positive values mean bright stuff in top channel, negative ones
    % mean bright stuff in bottom channel, and values around zero mean bright
    % stuff in both channels, ie., colocalised.
    % colocalis_threshold = 0.1; % values between -+colocalis_threshold are considered colocalised.
    % Plot histogram:
    subplot(3,2,2);
    hist(colocalis_values,100) % use 100 bins.
    xlabel('colocalisation value')
    ylabel('frequency')
    xlim([-0.8 0.8])
    
    % Colocalised fraction:
    % Percentage fraction of all values which are between -+colocalis_threshold:
    colocalised_fraction = 100*length(colocalis_values(abs(colocalis_values)<=colocalis_threshold))/length(colocalis_values);
    topChannel_fraction = 100*length(colocalis_values(colocalis_values>colocalis_threshold))/length(colocalis_values);
    bottomChannel_fraction = 100*length(colocalis_values(colocalis_values<(-colocalis_threshold)))/length(colocalis_values);
    
    {colocalised_fraction,topChannel_fraction,bottomChannel_fraction}
    
    subplot(3,2,4);
    axis off
    str1(1) = {['Image ',image_label,', Frame ',num2str(n)]};
    str1(2) = {['Colocalised fraction (%)= ',num2str(colocalised_fraction)]};
    str1(3) = {['TopChannel fraction (%)= ',num2str(topChannel_fraction)]};
    str1(4) = {['BottomChannel fraction (%)= ',num2str(bottomChannel_fraction)]};
    str1(5) = {['Colocalisation Threshold = ',num2str(colocalis_threshold)]};
    text(0,0.5,str1)
    
    % Sum of top and bottom channels instead:
    subplot(3,2,6);
    result3 = full_signalMask.*(Btop2_singleFrame + Bbot_singleFrame) - 2.*full_bgndMask;
    imshow(result3,[0,max(max(result3))]);
    colorbar('location','eastoutside'); title('shifted top plus bottom');
    colormap(jet);    
    
    % Save current figure (colour plot) as .png in output folder and close image:
    figName = strcat('plot_coloc',image_label,'_',bf_image_label,'_Frame',num2str(n),output_label); % name of figure file to save to.
    saveFigurePNG(output_folder_name,figName)
    
    % Save to output result structure;
    % save inputs first (same for all frames, but save anyway):
    result_struct(i).image_label = image_label;
    result_struct(i).bf_image_label = bf_image_label;
    result_struct(i).start_frame = start_frame;
    result_struct(i).end_frame = end_frame;
    result_struct(i).start_frame_forFrameAvg = start_frame_forFrameAvg;
    result_struct(i).end_frame_forFrameAvg = end_frame_forFrameAvg;
    result_struct(i).roi = roi;
    result_struct(i).output_label = output_label;
    result_struct(i).translateBy = translateBy; % amount top channel is shifted by with respect to bottom channel.
    % save results for current frame "n":
    result_struct(i).frameNumber = n;
    result_struct(i).colocalis_values = colocalis_values; % column vector.
    result_struct(i).colocalis_threshold = colocalis_threshold;
    result_struct(i).colocalised_fraction = colocalised_fraction;
    result_struct(i).topChannel_fraction = topChannel_fraction;
    result_struct(i).bottomChannel_fraction = bottomChannel_fraction;
    
    i = i+1; % advance saving index
end

[freqCounts binCentres]=hist(all_colocal_values,100);
scrsz = get(0,'ScreenSize'); % scrsz is a vector [left, bottom, width, height], left and bottom are equal to 1.
figure('Position',[scrsz(3)/3, scrsz(4)/3, scrsz(3)/1.5, scrsz(4)/2.5])
subplot(1,2,1)
bar(binCentres,freqCounts,'FaceColor',[0.8 0.8 0.8])
xlabel('colocalisation value (all frames)')
ylabel('frequency')
xlim([-1 1])
title(['mean ' num2str(mean(all_colocal_values)) '  stdev ' num2str(std(all_colocal_values))])
hold on;


% -----------------------
% % Fit peak to a Gaussian:

% fun_Gauss = fittype('A*exp(-(x-x0)^2/(2*sigma^2))','independent','x'); % define Gaussian funtion to fit to, with 'x' as independent variable.
% options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% % Guesses for fit parameters:
% guess_A = 350; % set by eye from plots.
% guess_sigma = 0.2; % set by eye from plots.
% guess_x0 = 0.2; % set by eye from plots.
% % Use coeffnames(fit_result) later to find out order of parameters.
% options.StartPoint = [guess_A guess_sigma guess_x0]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
% options.Lower = [ ]; % Lower bounds for fit parameters. In order: A, sigma, x0.
% options.Upper = [ ]; % Upper bounds for fit parameters. In order: A, sigma, x0.
% [fit_result gof] = fit(binCentres',freqCounts',fun_Gauss,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% % fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
% fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
% A_fit = fit_param_values(1); 
% sigma_fit = fit_param_values(2);
% x0_fit = fit_param_values(3);
% rsq_fit = gof.rsquare; % rsquare coefficient of fit.
% errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
% errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
% stDev_A = errorSTDEV(1);
% stDev_sigma = errorSTDEV(2);
% stDev_x0 = errorSTDEV(3);
% 
% disp('Gaussian fit to distribution:')
% disp(['A = ' num2str(A_fit) ' +- ' num2str(stDev_A)])
% disp(['x0 = ' num2str(x0_fit) ' +- ' num2str(stDev_x0)])
% disp(['sigma = ' num2str(sigma_fit) ' +- ' num2str(stDev_sigma)])
% disp(['rsq_fit = ' num2str(rsq_fit)])
% 
% % Plot fit:
% plot(fit_result,'k')
% hold off;


% -----------------------
% Fit peak to a sum of Gaussians:
% Sum of three Gaussians, one to the left (negative values (stuff on bottom
% channel)) centred around xleft and of width sigma_left; one to the right
% similarly, and one centred around zero, of width sigma_centre.

fun_Gaussians = fittype('A_left*exp(-(x-xleft)^2/(2*sigma_left^2))+A_centre*exp(-x^2/(2*sigma_centre^2))+A_right*exp(-(x-xright)^2/(2*sigma_right^2))','independent','x'); % define Gaussian funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_xleft = -0.2;
guess_A_left = freqCounts(abs(binCentres-(guess_xleft))==min(abs(binCentres-(guess_xleft))));  % 350
guess_sigma_left = 0.2;
guess_A_centre = freqCounts(abs(binCentres)==min(abs(binCentres)));  % 200
guess_sigma_centre = 0.2; 
guess_xright = 0.2; 
guess_A_right = freqCounts(abs(binCentres-(guess_xright))==min(abs(binCentres-(guess_xright))));  % 350
guess_sigma_right = 0.2;
% Order of parameters is: 'A_centre', 'A_left', 'A_right', 'sigma_centre',
% 'sigma_left', 'sigma_right', 'xleft', 'xright'.
% Use coeffnames(fit_result) later to find out order of parameters.
options.StartPoint = [guess_A_centre guess_A_left guess_A_right guess_sigma_centre guess_sigma_left guess_sigma_right guess_xleft guess_xright]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
% Constraints: A coeffs between 0 and Inf; sigma_centre between 0 and
% 0.5; sigma_left/right between 0.15 and 0.25; xleft between -0.4 and -0.1; xright beween 0.1 and 0.4. 
options.Lower = [0 0 0 0 0.15 0.15 -0.4 0.1]; % Lower bounds for fit parameters. In order.
options.Upper = [Inf Inf Inf 0.5 0.25 0.25 -0.1 0.4]; % Upper bounds for fit parameters. In order.
set(options,'Maxiter',6000); % 3000
set(options,'MaxFunEvals',8000); % 4000
[fit_result gof] = fit(binCentres',freqCounts',fun_Gaussians,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.

A_centre_fit = fit_param_values(1);
A_left_fit = fit_param_values(2); 
A_right_fit = fit_param_values(3); 
sigma_centre_fit = fit_param_values(4); 
sigma_left_fit = fit_param_values(5); 
sigma_right_fit = fit_param_values(6); 
xleft_fit = fit_param_values(7); 
xright_fit = fit_param_values(8); 

rsq_fit = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_A_centre_fit = errorSTDEV(1);
stDev_A_left_fit = errorSTDEV(2);
stDev_A_right_fit = errorSTDEV(3);
stDev_sigma_centre_fit = errorSTDEV(4);
stDev_sigma_left_fit = errorSTDEV(5);
stDev_sigma_right_fit = errorSTDEV(6);
stDev_xleft_fit = errorSTDEV(7);
stDev_xright_fit = errorSTDEV(8);

disp('----')
disp('Fit of distribution to sum of Gaussians:')
disp(['rsq_fit = ' num2str(rsq_fit)])
disp(['A_centre = ' num2str(A_centre_fit) ' +- ' num2str(stDev_A_centre_fit)])
disp(['sigma_centre = ' num2str(sigma_centre_fit) ' +- ' num2str(stDev_sigma_centre_fit)])
disp('----')
disp(['xleft = ' num2str(xleft_fit) ' +- ' num2str(stDev_xleft_fit)])
disp(['A_left = ' num2str(A_left_fit) ' +- ' num2str(stDev_A_left_fit)])
disp(['sigma_left = ' num2str(sigma_left_fit) ' +- ' num2str(stDev_sigma_left_fit)])
disp('----')
disp(['xright = ' num2str(xright_fit) ' +- ' num2str(stDev_xright_fit)])
disp(['A_right = ' num2str(A_right_fit) ' +- ' num2str(stDev_A_right_fit)])
disp(['sigma_right = ' num2str(sigma_right_fit) ' +- ' num2str(stDev_sigma_right_fit)])
disp('----')

% Plot fit:
% plot(fit_result,'k') % plot the sum of three Gaussians.

% Plot the three Gaussians separatedly too:
vector_x = binCentres'; % column vector.
step_x = vector_x(2)-vector_x(1);
vector_y_left = A_left_fit*exp(-(vector_x-xleft_fit).^2/(2*sigma_left_fit^2));
vector_y_centre = A_centre_fit*exp(-vector_x.^2/(2*sigma_centre_fit^2));
vector_y_right = A_right_fit*exp(-(vector_x-xright_fit).^2/(2*sigma_right_fit^2));

plot(vector_x,vector_y_left+vector_y_centre+vector_y_right,'k','LineWidth',1.5) % plot the sum of three Gaussians.

plot(vector_x,vector_y_left,'g','LineWidth',2) %  green, left, negative.
plot(vector_x,vector_y_centre,'b','LineWidth',2) % blue, central, zero.
plot(vector_x,vector_y_right,'r','LineWidth',2) % right, red, positive.
legend('','sum3Gaussians','green channel','both channels','red channel','Location','Best')

hold off;

% Calculate percentage fractions of area under each curve (left, right,
% central):
area_total = sum(step_x*(vector_y_left+vector_y_centre+vector_y_right));
area_left = sum(step_x*vector_y_left);
area_centre = sum(step_x*vector_y_centre);
area_right = sum(step_x*vector_y_right);
% Percentage fraction of total area corresponding to each Gaussian:
percent_frac_left = 100*area_left/area_total;
percent_frac_centre = 100*area_centre/area_total;
percent_frac_right = 100*area_right/area_total;
% Display on command window:
disp(['left fraction (%) = ',num2str(percent_frac_left)]);
disp(['centre fraction (%) = ',num2str(percent_frac_centre)]);
disp(['right fraction (%) = ',num2str(percent_frac_right)]);

subplot(1,2,2)
axis off
str1(1) = {['Image sequence ',image_label]};
str1(2) = {['Frames ',num2str(start_frame),' to ',num2str(end_frame)]};
str1(3) = {['Fit of distribution to sum of Gaussians: ','rsq-fit = ',num2str(rsq_fit)]};
str1(4) = {['A-centre = ',num2str(A_centre_fit),' +- ',num2str(stDev_A_centre_fit)]};
str1(5) = {['sigma-centre = ',num2str(sigma_centre_fit),' +- ',num2str(stDev_sigma_centre_fit)]};
str1(6) = {['xleft = ' num2str(xleft_fit) ' +- ' num2str(stDev_xleft_fit)]};
str1(7) = {['A-left = ' num2str(A_left_fit) ' +- ' num2str(stDev_A_left_fit)]};
str1(8) = {['sigma-left = ' num2str(sigma_left_fit) ' +- ' num2str(stDev_sigma_left_fit)]};
str1(9) = {['xright = ' num2str(xright_fit) ' +- ' num2str(stDev_xright_fit)]};
str1(10) = {['A-right = ' num2str(A_right_fit) ' +- ' num2str(stDev_A_right_fit)]};
str1(11) = {['sigma-right = ' num2str(sigma_right_fit) ' +- ' num2str(stDev_sigma_right_fit)]};
str1(12) = {['left fraction (%) = ',num2str(percent_frac_left)]};
str1(13) = {['centre fraction (%) = ',num2str(percent_frac_centre)]};
str1(14) = {['right fraction (%) = ',num2str(percent_frac_right)]};
text(0,0.5,str1)

% ----------------------------
% Save current figure (colour plot) as .png in output folder and close image:
figName = strcat('histogram',image_label,'_',bf_image_label,'_allvalues_Frames',num2str(start_frame),'to',num2str(end_frame),output_label); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName)

