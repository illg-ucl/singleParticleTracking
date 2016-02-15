function quickGraphColocalisation(image_label,bf_image_label,start_frame,end_frame,roi,output_label)
%
% Created by Isabel Llorente-Garcia, August 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Colocalisation analysis for frame averages.
%
% Produce a quick qualitative plot of possible colocalisation or
% antilocalisation of bright spots in the top and bottom channels of
% fluorescence images.
% See
% C:\Isabel\DataAnalysis\myAnalysis\quickGraphColocalisation\quickGraphColocalisation.doc
%
% Inputs:
% - image_label (fluorescence image sequence) such as '513', '490', etc... which corresponds to a certain .sif number of image sequence
% - bf_image_label: bright field image sequence (usually 10 frames only).
% - start_frame,end_frame: frames to use for the frame average. 
% "end_frame" can be set to 'end' when all frames in the image sequence, from start_frame to
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
%
% Example of how to call this function:
% quickGraphColocalisation('498','497',5,500,'full_image','A') for fluorescence image sequence 498, finding the displacement between top and bottom channels from bright field image '497'. .
% Another example: quickGraphColocalisation('498','497',5,'end',[10 60 10 60],''). Using
% last frame as end_frame.
%
% Note: it works also on a frame by frame basis, i.e., with start_frame =
% end_frame, though the thresholding for the signal mask might fail for one
% frame only if fluorescence is too dim...


%% Parameters: 

% CHECK:
% Colocalisation threshold:
% Positive final values mean bright stuff in top channel, negative ones
% mean bright stuff in bottom channel, and values around zero mean bright
% stuff in both channels, ie., colocalised.
colocalis_threshold = 0.144; % values between -+colocalis_threshold are considered colocalised.

% Parameter for thresholding the image to get the cell-signal mask within
% function getCellMaskAndBoundary3.m. The larger it its, the higher the
% threshold and the smaller the cell-signal mask region.
% If equal to 1, function getCellMaskAndBoundary3.m is equivalent to
% function getCellMaskAndBoundary.m
factor_forThresholding = 1.5; % (1.5 for oxphos data Isabel.)


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
[y1 x1] = overlapBF(bf_image_label,1,'end',30,30,30,30,0,0); % excluding 30 pixels on each edge.
translateBy = [y1 x1];
% ---------------------
% % Or choose by hand the amount by which to shift the image (down and left):
% translateBy = [6,-1];
% translateBy = [3,-5];
% translateBy = [2,-5];
% % ---------------------
% CHECK:
% % If you want to skip the translation:
% translateBy = [0 0];
% % ---------------------
disp('Top channel shifted by :');
disp(translateBy)


%% Calculate frame average between start_frame and end_frame:

disp('Fluorescence image:')
fr_avg = frameAverage(image_label,start_frame,end_frame,0,0); 


%% Separate top and bottom channels (usually red and green channels):

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


%%  Shift the top image:

% Shift the top channel to counteract slight displacement between channels:
% eg. [5,-4]: shift down 5 pixels and left 4 pixels.
% Shifted top image:
Btop02 = circshift(Btop0,translateBy); % shift top image by vector translateBy.
% "0" at end means full image before taking the roi.
result1 = Btop02-Bbot0; % shifted top image minus bottom one.
% The interesting plot is that of Btop02-Bbot0, with values between -1 and 1.
% Bright areas with values close to 1 will correspond to spots found in the
% top frame but not in the bottom one. 
% Dark areas with values close to -1 will correspond to spots found in the
% bottom frame but not in the top one.
% Grey areas with values close to 0 will correspond to either the
% background or areas where bright spots are colocalised in both channels.


%% Plot result in black and white, no cell/background masks used yet, full image before roi:

figure
subplot(3,1,1); imshow(Btop0,[]); colorbar('location','eastoutside'); title('shifted top channel'); % top greyscale image (0 to 1)
subplot(3,1,2); imshow(Bbot0,[]); colorbar('location','eastoutside');  title('bottom channel'); % bottom greyscale image (0 to 1)
subplot(3,1,3); imshow(result1,[]); colorbar('location','eastoutside'); title('shifted top minus bottom'); % shifted top image minus bottom one, greyscale.

% Save current figure (black and white plot) as .png in output folder and close image:
figName1 = strcat('plot_',image_label,'_',bf_image_label,'_frames',num2str(start_frame),'to',num2str(end_frame),output_label); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName1)


%% Use a region of interest from here on:
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


%% Re-scale images to make them between 0 and 1:

% Note that the re-scaling is done after shifting and taking the roi.
Btop2 = mat2gray(Btop2); % convert to grayscale.
Bbot = mat2gray(Bbot); % convert to grayscale.


%% Get signal and background masks (cell foreground shapes):

[SignalMaskTop2 cellBoundaryTop2]= getCellMaskAndBoundary3(Btop2,factor_forThresholding); % for shifted top channel, after taking roi.
% 1.5 is thte factor for the thresholding, the higher, the larger the
% threshold.
[SignalMaskBot cellBoundaryBot]= getCellMaskAndBoundary3(Bbot,factor_forThresholding); % for bottom channel, after taking roi.
BgndMaskTop2 = ~SignalMaskTop2; %  background mask is the complementary of the signal mask, shifted top channel, after taking roi.
BgndMaskBot = ~SignalMaskBot; % background mask is the complementary of the signal mask, bottom channel, after taking roi.

full_signalMask = (SignalMaskTop2 | SignalMaskBot); % union of signal masks of top and bottom channels.
full_bgndMask = (BgndMaskTop2 & BgndMaskBot); % intersection of background masks of top and bottom channels.

% CHECK:
% % % Alternative, to test threshold with single-label strains: use signal mask
% % % of only top channel and not its union with the bottom mask:
% full_signalMask = SignalMaskTop2;
% full_bgndMask = BgndMaskTop2;


%% Colocalisation result:

% result2 is the final result, for a colour plot. The background is now
% lower than any values in the foreground (-2), values close to +1 (or
% actually max1) correspond to spots found in the top frame but not in the bottom one,
% values close to -1 (or actually min1) correspond to spots found in the 
% bottom frame but not the top one, and values close to 0 correspond to
% background or co-localised spots.
% Note: (BgndMaskTop2 & BgndMaskBot) is the final background mask: intersection of the translated top one and the bottom one (ones in
% bgnd and zeros in cell/foreground region). 
% (SignalMaskTop2 | SignalMaskBot) is the union of the two signal cell
% masks.

result2_prelim = Btop2-Bbot; % preliminary result.
max1 = max(max(result2_prelim));
min1 = min(min(result2_prelim));

result2 = full_signalMask.*(Btop2-Bbot) - 2.*full_bgndMask;
result2_topPlot = full_signalMask.*Btop2 - full_bgndMask;
result2_botPlot = full_signalMask.*Bbot - full_bgndMask;

figure;
subplot(3,2,1); imshow(result2_topPlot,[0,1]); colorbar('location','eastoutside'); title('shifted top channel');
subplot(3,2,3); imshow(result2_botPlot,[0,1]); colorbar('location','eastoutside'); title('bottom channel');
subplot(3,2,5); imshow(result2,[-0.7,0.7]); % plot using a fixed colour/value range.
% subplot(3,2,5); imshow(result2,[-max(abs(min1),max1),max(abs(min1),max1)]); 
colorbar('location','eastoutside'); title('shifted top minus bottom');
colormap(jet); 


  

%% Output, colocalisation result:

colocalis_values = result2(full_signalMask==1); % Select only values inside mask. Colocalisation values (between -1 and 1).
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

subplot(3,2,4);
axis off
str1(1) = {['Image sequence ',image_label]};
str1(2) = {['Average of frames ',num2str(start_frame),' to ',num2str(end_frame)]};
str1(3) = {['Colocalised fraction (%)= ',num2str(colocalised_fraction)]};
str1(4) = {['TopChannel fraction (%)= ',num2str(topChannel_fraction)]};
str1(5) = {['BottomChannel fraction (%)= ',num2str(bottomChannel_fraction)]};
str1(6) = {['Colocalisation Threshold = ',num2str(colocalis_threshold)]};
text(0,0.5,str1)
% Save current figure (colour plot) as .png in output folder and close image:
figName2 = strcat('plotsColoc',image_label,'_',bf_image_label,'_frames',num2str(start_frame),'to',num2str(end_frame),output_label); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName2)

{colocalised_fraction,topChannel_fraction,bottomChannel_fraction}

% -------------------------------------
% Export colocalisation results:
cd(output_folder_name);
output_filename = strcat('colocValuesFrameAvg_',image_label,output_label,'.xls'); % output .xls filename. 
data_to_export = colocalis_values; % column vector.
warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.
xlswrite(output_filename,data_to_export,'coloc_values_onFrameAvg'); % write data to sheet 'coloc_values_onFrameAvg' in excel file.
% save vector of colocalisation values in a mat file too:
output_filename2 = strcat('colocValuesFrameAvg_',image_label,output_label,'.mat'); % output .mat filename.
save(output_filename2,'colocalis_values') % save variable colocalis_values within a .mat file named as output_filename2.

cd('..');
% ---------------------------------------

