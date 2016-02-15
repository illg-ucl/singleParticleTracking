function A_image_sequence = simulateImageSequence(image_label,frame_size,numFrames,photobleach_tau)
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

gaussian_sigma = 20; % gaussian width in pixels.
% keep a ratio of ~20 between frame_size/2 and gaussian_sigma.

% Noise parameters:
offset_noise_factor = 0.1; % For constant offset noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude (it is 1) to noise amplitude). 
shot_noise_factor = 0.5; % For noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity.


%% Generate simulated image sequence:

% initialise final image-sequence matrix to store things in:
A_image_sequence = zeros(frame_size,frame_size,numFrames);

% Loop through frames:
for i=1:numFrames
    
    % First generate a matrix with one Gaussian bright spot in it (top left quarter of final frame):
    Atop00 = fspecial('gaussian',frame_size/2,gaussian_sigma); % square matrix of size frame_size/2, with a normalised Gaussian spot of sigma "gaussian_sigma" pixels. This is the top left quarter of the final image.
    Atop0 = mat2gray(Atop00); % scale between 0 and 1 (to grayscale image).
    
    % Account for exponential photobleaching decay;
    bleach_factor = exp(-(i-1)/photobleach_tau);
    Atop = bleach_factor*Atop0;
    
    % imshow(Atop,[])
    % plot(1:frame_size/2,Atop(frame_size/4,:)) % show linear Gaussian
    % profile.
       
    % Generate bottom half of image:
    Abottom = zeros(frame_size/2,frame_size); % rectangular matrix, bottom half of final image.
    
    % Put pieces together into final frame image (noiseless so far):
    A = [[Atop zeros(frame_size/2,frame_size/2)]; Abottom];
    
    % Create noise for images:
    % Shot-noise: proportional to square root of number of photons/intensity:
    shot_noise_matrix = shot_noise_factor*sqrt(A).*rand(size(A,1),size(A,2));
    % Constant offset noise:
    constant_offset_noise_matrix = offset_noise_factor*rand(size(A,1),size(A,2));
    
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