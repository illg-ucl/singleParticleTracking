function A_image_sequence = simulateImageSequence(image_label,frame_size,numFrames,photobleach_tau_spot)
%
% Created by Isabel Llorente Garcia. April 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Create simulated image sequence (.mat file).
% Create image sequence of stacked 2D images and save as a 3D matrix in a mat file.
% In particular, this function creates an image with a bright spot on top
% of a background, and adds both offset noise and intensity-dependent noise
% to it. Only one bright spot in the top left quarter of the images is created.
%
% INPUTS:
% - image_label: label to save the image as image_sequence_ followed by the
% label. The label can be for example '1', so that image is saved as image_sequence_1. 
% - frame_size: size of final frames in sequence (squared ones). 
% - numFrames: number of frames in final simulated image sequence.
% - photobleach_tau: exponential photobleaching decay time constant, in
% frames (not seconds!). To simulate photobleaching from frame to frame.
% For time between frames of 40ms, a time constant tau of 1 second
% corresponds to 25 frames, tau = 2s corresponds to 50 frames, etc.
% Default: 25-75 frames.
%
% OUTPUTS:
% A .mat file is saved with a MxMxJ matrix, where MxM is the size of each
% frame (M=frame_size), and J=numFrames is the number of frames in the image sequence.


%% PARAMETERS:

gaussian_sigma = 3; % gaussian width in pixels.
autoflu_radius = 50; % radius of the circle mimicking the cell and its autofluorescence n.b. 1 extra pixel is added to total diameter by the fspecial function used
% keep a ratio of ~20 between frame_size/2 and gaussian_sigma.

% Noise parameters:
offset_noise_factor = 0.1; % For constant offset noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude to noise amplitude). 
shot_noise_factor = 0.5; % For noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity.
cell_autofluo_noise_stoch_factor = 0.5; % - Stochastic element of noise mimicking autofluorescence of the cell  %G 
cell_autofluo_intensity = 0.15; % - mean intensity of  autofluorescence in the simulated cell before photobleaching (scale 0-1) %G
photobleach_tau_autofluo = 18; % rate of bleaching of cell autofluorescence

%% Generate simulated image sequence:

% initialise final image-sequence matrix to store things in:
A_image_sequence = zeros(frame_size,frame_size,numFrames);

    
% First generate a matrix with one Gaussian bright spot in it (top left quarter of final frame):
Atop00 = fspecial('gaussian',frame_size/2,gaussian_sigma); % square matrix of size frame_size/2, with a normalised Gaussian spot of sigma "gaussian_sigma" pixels. This is the top left quarter of the final image.
Atop0 = mat2gray(Atop00); % scale between 0 and 1 (to grayscale image).

% Loop through frames:
for i=1:numFrames
    
    % Account for exponential photobleaching decay (of spot);
    bleach_factor = exp(-(i-1)/photobleach_tau_spot);
    Atop = bleach_factor*Atop0;
    
    % imshow(Atop,[])
    % plot(1:frame_size/2,Atop(frame_size/4,:)) % show linear Gaussian
    % profile.
    if i==1
        max_ampli = max(max(Atop)); % maximum value, gaussian amplitude.
    end
    
    % Generate bottom half of image:
    
    %Create uniform circle to simulate the autofluorescence from the cell
    
    circle = fspecial('disk', autoflu_radius); % Creates uniform intensity circle of radius 'autoflu_radius' on a black square background
    scaled_circle = mat2gray(circle)*cell_autofluo_intensity; % scales between 0 and 1 (to greyscale image) and then scales for the intensity parameter
    total_autoflu_diameter = size(scaled_circle,1); % gets the overall size of the square containing the uniform circle
    
    % add andom noise in autofluorescence (bigger than a pixel)
    num_of_noise_spots=40; % this is the number of bright, multi-pixel autofluorescence spots per area of (total_autoflu_diameter)squared
    noise_spots_matrix=zeros(total_autoflu_diameter); % creates area in which to place these spots
    for j=1:num_of_noise_spots % creates random coordinates for the spots
        x = randi(total_autoflu_diameter);
        y = randi(total_autoflu_diameter);
        noise_spots_matrix(x,y)=1;
    end
    filter_1=fspecial('disk', 2); % filter to enlarge spot around assigned coordinates
    spots_with_filter_1=imfilter(noise_spots_matrix,filter_1); % apply this filter
    scaled_spots_with_filter_1=mat2gray(spots_with_filter_1); % scale the matrix produced by this filter between 0 and 1
    filter_2 = fspecial('gaussian',4,3); % filter to blur the edges of these spots
    spots_with_filter_2 = imfilter(scaled_spots_with_filter_1,filter_2); %apply this filter
    spots_in_circle = 0.8*spots_with_filter_2.*scaled_circle; % exclude all spots assigned to coordinates outside the cell 
    scaled_circle = scaled_circle + spots_in_circle; % add to the uniform autofluorescence
    
    % add random noise in autofluorescence (at single pixel level)
    Abottom00 = scaled_circle+scaled_circle.*(cell_autofluo_noise_stoch_factor*(-1+rand(total_autoflu_diameter, total_autoflu_diameter)*2)); % adds/subtracts normally distributed random amount of autofluorescence
    
    % makes a black 'border' around autofluorescence disk, so that the square has an area  one quarter the size of the final image
    width_padding_pre_matrix = (frame_size/2-total_autoflu_diameter)/2 + 0.5; % this line and the next line finds the width of the 'border' needed to be concatenated around the square matrix containing the circle. The +/- 0.5 is so that the border width is an integer number of pixels (as total_autoflu_diameter, which is always odd, has been divided by 2) 
    width_padding_post_matrix = (frame_size/2-total_autoflu_diameter)/2 - 0.5; % see above line
    padded_pre_matrix = padarray(Abottom00,[width_padding_pre_matrix width_padding_pre_matrix], 'pre'); % adds the 'border' at the top and left
    padded_post_matrix = padarray(padded_pre_matrix,[width_padding_post_matrix width_padding_post_matrix], 'post'); % adds the 'border' at the bottom and right
    
    % Account for exponential photobleaching decay of cell autofluorescence and put bottom half of image together;
    bleach_factor = exp(-(i-1)/photobleach_tau_autofluo);
    Abottom0 = bleach_factor*padded_post_matrix;
    Abottom = horzcat(Abottom0, zeros(frame_size/2,frame_size/2)); % puts together bottom half of the image, so the autofluorescence circle is in the bottom left quarter 
    
    % Overlay the cell autofluorescence and the Gaussian bright spot for top image
    Atop_with_autofluo = Abottom0 + Atop; 
    
    % Put pieces together into final frame image (noiseless except for the autofluorescence specific noise):
    A = [[Atop_with_autofluo zeros(frame_size/2,frame_size/2)]; Abottom];
    A = max(min(A,1),0);  %puts upper limit on any intensity values, saturation of pixels at value 1
    
    % Create noise for images:
    % Shot-noise: proportional to square root of number of photons/intensity:
    shot_noise_matrix = shot_noise_factor*sqrt(A).*rand(size(A,1),size(A,2));
    % Constant offset noise:
    constant_offset_noise_matrix = max_ampli*offset_noise_factor*rand(size(A,1),size(A,2));
    
    % Include noise on images:
    A_with_noises = A + shot_noise_matrix + constant_offset_noise_matrix;
    
    % imshow(A_with_noises,[])
    % figure;
    % plot(1:2*matrix_size,A_with_noises(matrix_size/2,:)) % see linear
    % Gaussian profile.
    
    % stack the generated frames onto a frame_size-by-frame_size-by-numFrames
    % matrix, which is the image sequence:
    A_image_sequence(:,:,i) = A_with_noises;
    
end

%% Save final image sequence:

% Save image sequence as a matrix in a .mat file with the corresponding image label:
save(strcat('image_sequence_',image_label),'A_image_sequence')


% Graphic test, display image sequence as a video at the end:
for i=1:numFrames
imshow(A_image_sequence(:,:,i),[0,1.5]); 
pause(0.2);
end
