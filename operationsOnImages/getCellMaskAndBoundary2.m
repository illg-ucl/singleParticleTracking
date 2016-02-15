function [SignalMask CellBoundaryMask] = getCellMaskAndBoundary2(frame,local_region)
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Same as function getCellMaskAndBoundary.m but using a local region to find the thresholding value.
% And not dividing image in two parts.
% Use correlation of an image with a small gaussian mask and gradient to
% obtain a mask with ones at cell positions and zeros in background
% (SignalMask, binary image), and the CellBoundaryMask which is another
% binary image with ones at cell boundaries and zeroes everywhere else.
%
% Inputs: "frame" is the input image. It must have values between 0 and 1.
% IMPORTANT NOTE: Input image frame must have values between 0 and 1. 
% Use mat2gray(frame) instead of frame if values are not between 0 and 1.
% IMPORTANT NOTE: input image is a full image.
%
% local_region: region to use to find threshold of image. It is given as
% [xleft xright ytop ybottom] with all values in pixels.
%
% Note: functions "dftcorr" and "imgrad" can be found in
% C:\Isabel\myMatlabFiles.
wcorr = fspecial('gaussian',10,5); % builds Gaussian mask of size 10 by 10 and sigma equal to 5 pixels.
% Note that it is better if the size of the wcorr matrix is an even number (for later).
gcorr = dftcorr(frame,wcorr); % dftcorr calculates the correlation of the image with the mask 
% (it finds places in the image frame which match the mask wcorr).

% Local region to use for thresholding the image: 
% local_region = [xleft xright ytop ybottom].
xleft = local_region(1);
xright = local_region(2);
ytop = local_region(3);
ybottom = local_region(4);

my_threshold = graythresh(gcorr(ytop:ybottom,xleft:xright));
% my_threshold = graythresh(gcorr(half+1:2*half,:)); % This is how it was
% in getCellMaskAndBoundary.m.

% To separate foreground and background we threshold the correlation image:
sigBW = im2bw(gcorr,my_threshold); % thresholding 
% Join top and bottom thresholded images into a single image of the original size:
gBW = sigBW; % thresholded correlation image.
gBW2 = circshift(gBW,[round(size(wcorr,1)/2) round(size(wcorr,2)/2)]); % This shifts the image down and to the right by half the size of the mask wcorr.
% The shifting is necessary because the correlation image gives the
% positions of the elements that match the mask as the top left corner
% position of the mask.

% This is the first function output:
SignalMask = gBW2; % This is the foreground (or signal) mask with ones inside cells and zeros in background. It is a binary image.

% Now obtain cell boundaries:
SignalMask2 = mat2gray(SignalMask); % First need to convert binary SignalMask to a grayscale image to be able to apply gradient.
cellBoundaries = imgrad(SignalMask2); % Grayscale image with cell boundaries.
cellBoundaries2 = im2bw(cellBoundaries); % Convert cellBoundaries to a binary image with ones at cell edges or boundaries.
CellBoundaryMask = cellBoundaries2; % This is the second function output.


%---------------------------------------------------------------------- 

% % Check result graphically: overlay found boundaries on original image frame:
% subplot(2,3,1); 
% imshow(frame,[]); 
% title('Original image');
% 
% subplot(2,3,2); 
% imshow(SignalMask,[]); 
% title('Signal Mask');
% 
% subplot(2,3,3); 
% imshow(CellBoundaryMask,[]); 
% title('cell-boundaries Mask');
% 
% subplot(2,3,4); 
% imshow(frame.*SignalMask,[]); 
% title('Signal Mask * original image');
% 
% subplot(2,3,5); 
% imshow(frame.*(~SignalMask),[]); % the background mask would be the complement of the Signal mask.
% title('Background Mask * original image');
% 
% subplot(2,3,6); 
% imshow(CellBoundaryMask+frame,[]); 
% title('cell boundaries overlaid on original image');
% 
% 
