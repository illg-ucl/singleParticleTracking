function [SignalMask, Boundary] = getBrightSpotsMask(frame,factor_forThresholding)
%
% Created by Isabel Llorente-Garcia, August 2014.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Based on getCellMaskAndBoundary.m, with the possibility of adjusting
% the thresholding with a factor "factor_forThresholding".
% Use correlation of an image with a small gaussian mask and gradient to
% obtain a mask with ones at bright area positions and zeros in background
% (SignalMask, binary image), and the Boundary which is another
% binary image with ones at boundary (edge of bright region) and
% zeroes everywhere else.
%
% Inputs: "frame" is the input image. It must have values between 0 and 1.
% "factor_forThresholding": the larger it its, the higher the
% threshold and the smaller the cell-signal mask region.
% IMPORTANT NOTE: Input image frame must have values between 0 and 1. 
% Use mat2gray(frame) instead of frame if values are not between 0 and 1.
% IMPORTANT NOTE: no splitting in top and bottom channels, same operations
% are carried out for the entire image.
%
% Note: functions "dftcorr" and "imgrad" can be found in
% C:\Isabel\myMatlabFiles.
wcorr = fspecial('gaussian',10,5); % builds Gaussian mask of size 10 by 10 and sigma equal to 5 pixels.
% Note that it is better if the size of the wcorr matrix is an even number (for later).
gcorr = dftcorr(frame,wcorr); % dftcorr calculates the correlation of the image with the mask 
% (it finds places in the image frame which match the mask wcorr).

% To separate foreground and background we threshold the correlation image:
% Threshold entire image in one go:
gBW = im2bw(gcorr,factor_forThresholding*graythresh(gcorr)); % thresholding with adjusting factor 
gBW2 = circshift(gBW,[round(size(wcorr,1)/2) round(size(wcorr,2)/2)]); % This shifts the image down and to the right by half the size of the mask wcorr.
% The shifting is necessary because the correlation image gives the
% positions of the elements that match the mask as the top left corner
% position of the mask.

% This is the first function output:
SignalMask = gBW2; % This is the foreground or signal mask with ones inside bright area and zeros in background. It is a binary image.

% Now obtain boundary:
SignalMask2 = mat2gray(SignalMask); % First need to convert binary SignalMask to a grayscale image to be able to apply gradient.
Boundary_0 = imgrad(SignalMask2); % Grayscale image with bright region boundaries.
Boundary = im2bw(Boundary_0); % Convert Boundary image to a binary image with ones at bright-region edges or boundaries.
% Boundary is the second function output.


%---------------------------------------------------------------------- 

% Check result graphically: overlay found boundaries on original image frame:
subplot(2,3,1); 
imshow(frame,[]); 
title('Original image');

subplot(2,3,2); 
imshow(SignalMask,[]); 
title('Signal Mask');

subplot(2,3,3); 
imshow(Boundary,[]); 
title('bright-region Boundary');

subplot(2,3,4); 
imshow(frame.*SignalMask,[]); 
title('Signal Mask * original image');

subplot(2,3,5); 
imshow(frame.*(~SignalMask),[]); % the background mask would be the complement of the Signal mask.
title('Background Mask * original image');

subplot(2,3,6); 
imshow(Boundary,[]); 
title('Boundary overlaid on original image');


