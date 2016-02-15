function spot = calculateSpotIntensity(frame_data,x0,y0,subarray_halfwidth,inner_circle_radius)
% 
% Created by Isabel Llorente-Garcia, October 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Function to calculate intensity at a fixed position within a given frame,
% using a square subarray and a circular mask centred on fixed position
% (x0,y0). (So no tracking or finding spot centre involved.)
% 
% 
% INPUTS:  
% - frame_data is a matrix containing the image sequence data for a single frame. 
% frame_data is a matrix with original values, so not the grayscale version between 0 and 1.
% frame_data is later made of class double.
% It can be the output of extract1frame(frame_number) or
% extract_image_sequence_data(image_label), for example.
%
% - x0, y0 are fixed positions, centre of the bright spot. 
%
% - subarray_halfwidth: the algorithm selects this number of pixels above and below and to
% left and right of estimated centroid pixel to form the selected image
% subarray (I). (The default value is 8 pixels).
%
% - inner_circle_radius: radius of inner circular mask in pixels (default is 5 pixels). This mask moves
% around inside the subarray to find the bright spot.
%
% OUTPUTS: values resulting from the iterative method. The output is a
% structure "spot" with the following fields (see end of file).
%
% Example of how to call this function: s1 = calculateSpotIntensity(A,271,364,8,5), 
% where A is the image array (single frame data). Then to obtain results
% from the structure s1, do, for example: s1.CentreX, s1.CentreY, etc.

% make frame_data of class double:
frame_data = double(frame_data);

d = subarray_halfwidth; % subarray halfsize. d needs to be at least equal to the radius of the inner mask, inner_circle_radius. Default: d=8, inner_circle_radius=5.

%-----------------------
% error control:
if subarray_halfwidth < inner_circle_radius
    error('findSpotCentre1frame:one','check input parameters: "subarray_halfwidth" cannot be smaller than "inner_circle_radius"')
end

% If the spot candidate is at the edge of the image move it away from the edge until we can take a
% subarray around it:
if (round(y0)-d)<1
    y0 = d+1;
end
if (round(y0)+d)>size(frame_data,1)
    y0 = size(frame_data,1)-d;
end
if (round(x0)-d)<1
    x0 = d+1;
end
if (round(x0)+d)>size(frame_data,2)
    x0 = size(frame_data,2)-d;
end    
%------------------------

% Create image subarray (I) of size (2*d+1)x(2*d+1) (default size 17x17) centred on the
% centroid estimate pixel (x0,y0). I is fixed during the iterative process of finding the centre of the bright spot.
I = frame_data(round(y0)-d:round(y0)+d,round(x0)-d:round(x0)+d);

% Create matrices containing the x and y positions corresponding to the
% subarray I. Xs is a matrix of the same size as I, containing x values for
% all pixels and similarly for Ys.
[Xs,Ys] = meshgrid(round(x0)-d:round(x0)+d,round(y0)-d:round(y0)+d);


% Create inner circle mask: circle of radius inner_circle_radius pixels, with ones inside and zeros outside.
% This mask can move around inside the fixed image subarray I.
inner_mask = zeros(size(I,1)); % pre-define inner_mask size.
% radius of inner circle mask is inner-radius. Default is 5 pixels.
for ii = 1:size(I,1) % ii = index for rows
    for jj = 1:size(I,2) % jj = index for columns
        % inner circle mask: if distance to point (x0,y0) is <=inner_circle_radius, then 1, otherwise, 0:
        sum_sq = (Xs(ii,jj)-x0)^2 + (Ys(ii,jj)-y0)^2;
        if round(sqrt(sum_sq))<=inner_circle_radius
            inner_mask(ii,jj)=1;
        else
            inner_mask(ii,jj)=0;
        end
    end
end

% The background mask is just the negative of the inner circle mask, i.e.,
% not-inner_mask_0. It has zeros in centre and ones in surrounding pixels.
bgnd_mask = double(~inner_mask); % make of class double (not binary) for what follows.

pos_bgnd = find(bgnd_mask==1); % positions of bgnd only intensities in bgnd mask.
Ibgnd = I.*bgnd_mask; % bgnd region.

% ----------------------------------------
% Total background intensity in background region before background correction:
% Ibg_tot = sum(sum(I.*bgnd_mask)); % this gives same result.
Ibg_tot = sum(Ibgnd(pos_bgnd));
% Median background intensity in bgnd region, per pixel, Ibg_avg.
% The number of pixels in the background mask equal to 1 is
% length(find(bgnd_mask~=0)):
% Ibg_avg = median(Ibgnd(pos_bgnd)); % use median instead of mean for the bgnd to exclude hot pixels or nearby bright spots.
Ibg_avg = mean(Ibgnd(pos_bgnd)); % use mean for the bgnd to exclude hot pixels or nearby bright spots. This is the offset noise level per pixel before background subtraction.
% Total intensity in inner mask region before background correction:
Iinner_tot = sum(sum(I.*inner_mask));
%
% Calculate background-corrected subarray image:
I2 = I-Ibg_avg;

% Calculate standard deviation of remaining noise in background
% corrected image I2:
bg_noise_offset_afterBGsubtract = mean(I2(pos_bgnd)); % offset noise level per pixel after background subtraction, should be close to zero.
bg_noise_std = std(I2(pos_bgnd)); % standard deviation of matrix elements in bgnd region.
% Total spot intensity (within the inner circle mask), background corrected:
Isp = sum(sum(I2.*inner_mask));
% -----------------------------------------
   

% -----------------------------------------------
% Output of the function:
% The output is a structure "spot" with the following fields:
spot.CentreX = x0; % x-centre result found
spot.CentreY = y0; % y-centre result found
spot.IspTot = Isp; % Total intensity within the spot (inner circle) mask in subarray, background corrected.
spot.bg_noise_offset_afterBGsubtract = bg_noise_offset_afterBGsubtract; % offset background noise level after background subtraction, per pixel.
spot.BgNoiseStd = bg_noise_std; % standard deviation of matrix elements in bgnd region only, after background subtraction, per pixel.
spot.IbgAvg = Ibg_avg; % Average (per pixel) background intensity offset in original image subarray, before background correction.
spot.IbgTot = Ibg_tot; % Total background intensity in background region in original image subarray, before background correction.
spot.IinnerTot = Iinner_tot; % Total intensity inside inner mask region in original image subarray, before background subtraction.
spot.TrajNumber = []; % add field "trajectory number" for later, for function FindTrajects.m.
% -----------------------------------------------


% -----------------------------------------------------------------------
% Auxiliary stuff below:

% Plot things:
% figure
% subplot(2,3,1); imshow(I,[],'InitialMagnification',1000); title('original image subarray');
% subplot(2,3,2); imshow(inner_mask,[],'InitialMagnification',1000); title('inner circle mask');
% subplot(2,3,3); imshow(bgnd_mask,[],'InitialMagnification',1000); title('background mask');
% subplot(2,3,4); imshow(I2,[],'InitialMagnification',1000); title('bgnd corrected image');
% subplot(2,3,5); imshow(mask_gauss,[],'InitialMagnification',1000); title('gaussian mask');
% subplot(2,3,6); imshow(I2.*mask_gauss,[],'InitialMagnification',1000); title('gauss*bgnd-corrected');
%
% % Images of results:
% figure
% imshow(frame_data,[],'InitialMagnification',150);
% hold;
% plot(x0,y0,'o','Color','r','MarkerSize',8); title('obtained centre')
% hold off;
% 
% figure
% imshow(I,[],'InitialMagnification',1000);
% hold;
% plot(x0-x0+d+1,y0-y0+d+1,'o','Color','r','MarkerSize',70); title('obtained centre')
% hold off;