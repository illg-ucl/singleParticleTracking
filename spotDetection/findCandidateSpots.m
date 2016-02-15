function [xx yy] = findCandidateSpots(frame,method)
% 
% Created by Isabel Llorente-Garcia, August 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Find candidate (fluorescence) bright spots on images (of cells).
% The candidate spots can later be used as input of function
% "findSpotCentre1frame".
%
% Note: image thresholding occurrs in two halves: separating top and bottom
% halves.
%
% INPUTS:
% frame: Input image, a single frame from a recorded sequence. 
% IMPORTANT NOTE: Input image frame must have values between 0 and 1. 
% Use mat2gray(frame) instead of frame if values are not between 0 and 1.
%
% method: two different methods, either 1 or 2.
%
% OUTPUTS:
% xx: columns of positions of spot candidates
% yy: rows of positions of spot candidates
%
% Note that the input image can be obtained as:                       
% frame = extract1frame(framenumber); with framenumber=100, for example. 
%
% It can be useful to apply a Gaussian convolution to filter the image before
% using this function (e.g.  frame2 =
% imfilter(frame,fspecial('gaussian',3))).

frame0 = frame; % original image.

% To enable this part add input: preprocess: either 1 or 0. If equal to 1, image is preprocessed by
% calculating its gradient and opening it, and then methods 1 or 2 are applied.
% % Do some pre-procesing on the original image "frame":
% if preprocess==1
%     h=fspecial('sobel'); % filter mask to calculate gradient: h=[1 2 1; 0 0 0; -1 -2 -1] and h'=[1 0 -1; 2 0 -2; 1 0 -1].
%     % Calculate gradiente image:
%     grady=imfilter(frame,h,'replicate'); % vertical gradient.
%     gradx=imfilter(frame,h','replicate'); % horizontal gradient.
%     Grad=sqrt(grady.^2+gradx.^2); % gradient image.
%     Grad2=imopen(Grad,strel('disk',3)); % opened grad image to highlight bright spot regions.
%     frame = Grad2;
% end


% PARAMETERS:
disk_radius = 5; % disk radius in pixels for top hat transformation, quite critical to the method for finding candidate bright spots.
% ----------------------

se = strel('disk',disk_radius); % structural element, disk of radius 5 pixels (default). The radius is quite critical to the method.
Signal = imtophat(frame,se); % enhanced signal-image via "top hat" transformation (evens out the background). 
% "Top-hat": substract an opened image from the original. Opening with a
% large enough structural element (disk of radius at least 5 pixels) ensures that
% the opening of the image produces a reasonable estimate of the background
Bgnd0 = imopen(frame,se); % Background estimate. Morphological opening (=erosion followed by dilation) with large enough structural element ("se" does not fit entirely within the foreground).

% Middle row in image that separates top half (red channel) from bottom
% half (green channel). If the image is 512x512, half = 256.
half = round(size(frame,1)/2); % round to closest integer.

% Separating foreground (signal) from background:

% Thresholding each half:
SignalMaskTop = im2bw(Bgnd0(1:half,:),graythresh(Bgnd0(1:half,:))); % graythresh gives the threshold and im2bw converts image into a thresholded black and white (bw) image.
SignalMaskBottom = im2bw(Bgnd0(half+1:size(frame,1),:),graythresh(Bgnd0(half+1:size(frame,1),:))); % thresholding the bottom half.
% Join top and bottom masks into a single image of the original size:
SignalMask = [SignalMaskTop; SignalMaskBottom];

% % Background-only image:
% Bgnd = A.*imcomplement(SignalMask);
% imshow(Bgnd,[]);

% Signal (foreground) image:
% Signal2 = A.*SignalMask; % foreground or signal region.
Signal3 = Signal.*SignalMask; % enhanced signal-only image.

% imshow(SignalMask,[]);figure;

% Thresholding the signal-only image to try and find spot candidates:
sigBWtop = im2bw(Signal3(1:half,:),graythresh(Signal3(1:half,:))); % thresholding the top half
sigBWbottom = im2bw(Signal3(half+1:size(frame,1),:),graythresh(Signal3(half+1:size(frame,1),:))); % thresholding the bottom half
% Join top and bottom thresholded images into a single image of the original size:
sigBW = [sigBWtop; sigBWbottom];

B1 = imopen(sigBW,strel('disk',1)); % morph. opening with disk of radius 1 pixel.
B2 = bwmorph(B1,'fill',1); % fill holes of size one pixel, apply operation one time (1).
B3 = bwmorph(B2,'open',1); % open image, apply operation one time.


% Method 1:
if method==1
    B4 = imerode(B3, [1 1; 1 1]); % erode with a 2x2 square shape.
    B5 = imerode(B4, [1 1; 1 1]); % erode with a 2x2 square shape.
    B6 = bwmorph(B5,'shrink',1); % shrinks objects to points, apply operation one time.
    % B6 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
    result = B6;
end


% Method 2:
if method==2
    % B7 = bwdist(~B3); % distance transform (distance between each pixel and its nearest nonzero pixel) of the complement of the image
    B8 = bwulterode(B3); % "ultimate erosion" of a binary image = regional maxima of the Euclidean distance transform of the complement of the image.
    % subplot(2,2,1);imshow(frame,[]);
    % subplot(2,2,2);imshow(B3,[]);
    % subplot(2,2,3);imshow(B7,[]);
    % subplot(2,2,4);imshow(B8,[]);
    result = B8; % B8 is a binary image (matrix) with ones at positions of candidate spots and zeros elsewhere.
end

% Create matrices containing the x and y positions corresponding to the
% subarray I. Xpos is a matrix of the same size as I, containing x values for
% all pixels and similarly for Ypos.
[Xpos Ypos] = meshgrid(1:size(frame,2),1:size(frame,1));

% The positions of the candidate spots are:
xx=nonzeros(Xpos.*result); % "nonzeros" returns a column vector with the non-zero matrix elements ordered by columns.
yy=nonzeros(Ypos.*result);

% Graphically show the result overlaying result on original image:
% imshow(frame0+result,[]);
