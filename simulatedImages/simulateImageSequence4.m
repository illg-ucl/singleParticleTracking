function A_image_sequence = simulateImageSequence4(image_label,frame_size,numFrames,num_molecules,spot_distance)
%
% Started by Isabel Llorente-Garcia. August 2014.
%
% Create simulated image sequence (.mat file). Similar to
% simulateImageSequence3, but no bleaching, no padding, and it simulates two spots
% which are close together. Shot noise calculation uses poissrnd instead.
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
% - num_molecules % number of single-molecules (Gaussian intensity spots with amplitude I0 equal to 1) within one bright spot.
% - spot_distance: x-distance between bright spot centres, in pixels.
%
% OUTPUTS:
% A .mat file is saved with a MxMxJ matrix, where MxM is the size of each
% frame (M=frame_size), and J=numFrames is the number of frames in the image sequence.
%
% EXAMPLE of how to call this function: 
% A = simulateImageSequence3('try1',100,200,5);
%
% Note: gaussian bright spots have an amplitude of 1, noise added on
% top of that. 


%% PARAMETERS:

% For bright spots:
gaussian_sigma = 10; % Gaussian width of both spots in pixels (doesn't need to be an integer). (Default = 3.5).
num_bright_spots = 2; % Number of bright spots generated.
spot_position_i = 50; % x-value of position of first spot; also, y-value of position of first and second spots.
% The second spot is located at x = spot_position_i+spot_distance, and
% y=spot_position_i.

% Noise parameters:
% Note: all referred to the amplitude of a single-molecule bright spot Gaussian,
% which is 1. 
shot_noise_factor = 1; % Controls amplitude of shot noise (Poissonian noise), of noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity. (Default = 1).
offset_level = 0.2; % Offset background level. Default: 0.2.
bgnd_noise_std = 0.5; % For background noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude to noise standard deviation). (Default = 0.5). 


%% Error control:
% Make sure spots are not generated at positions that are outside the size
% of the frame (frame_size):
if round(frame_size) < (spot_position_i + spot_distance + 3*gaussian_sigma)
    disp ('----  !!!!!!  ----')
    disp('ERROR!! - Check parameters inside "simulateImageSequence4.m": check spot_position_i, spot_distance and gaussian_sigma relative to frame_size input value !!!')
    disp('Exiting function.')
    disp ('----  !!!!!!  ----')
    return
end

%% Generate simulated image sequence:

% Initialise final image-sequence matrix to store things in:
A_image_sequence = zeros(frame_size,frame_size,numFrames);

% GENERATE SIMULATED BRIGHT FLUORESCENT SPOTS:
% Generate two Gaussian spots in a fixed placed, horizontally separated by
% spot_distance;
% Generate a small matrix with one Gaussian bright spot in it (to act as filter later on, to simulate flu bright spots):
gauss_bright_spot = fspecial('gaussian',min(round(10*gaussian_sigma),frame_size),gaussian_sigma); % square matrix of size 10*gaussian_sigma, with a normalised Gaussian spot of sigma "gaussian_sigma" pixels. 
% Generate fixed positions for bright spots inside top-left quarter of image:
Xpos_bright_spots = zeros(num_bright_spots,1); % initialise column vector.
Ypos_bright_spots = zeros(num_bright_spots,1); % initialise column vector.
% Initialise top-left quarter of image with bright spots, for first frame:
A_spots_1 = zeros(frame_size,frame_size);
for k=1:num_bright_spots
    % Use fixed positions of choice:
    Xpos_bright_spots(k,1) = spot_position_i + spot_distance*(k-1);
    Ypos_bright_spots(k,1) = spot_position_i;

    A_spots_1(Ypos_bright_spots(k,1),Xpos_bright_spots(k,1)) = 1; % place a 1 at the positions of bright spots in A_top_left_quarter_1b.
end

disp('Bright spot positions x,y(pix):')
disp([Xpos_bright_spots Ypos_bright_spots])

% Place gaussian bright flu spots at positions given by above places, use filter:
A_spots_2 = imfilter(A_spots_1,gauss_bright_spot);
% Scale between 0 and 1 (to grayscale image), and multiply by the number of molecules per bright spot (each spot has a maximum Gaussian amplitude of 1):
A_spots_3 = num_molecules*mat2gray(A_spots_2); 
% Asume I0, the amplitude of a Gaussian intensity spot from a single
% molecule is 1. Then assume width stays the same when increasing the
% number of molecules and that I0 scales with the number of molecules,
% num_molecules, in one bright spot.


% LOOP THROUGH FRAMES:
for i=1:numFrames
     
    % Add noise to each frame on image sequence:
    % Shot-noise: proportional to square root of total number of photons/intensity:
    % The uncertainty in the no. of counts (standard deviation of the
    % noise) is proportional to the square root of the number of counts at
    % each position. Shot noise follows a Poissonian distribution:
    shot_noise_matrix = shot_noise_factor*poissrnd(A_spots_3);
    % Function poissrnd: A_spots_3 is the image with my Gaussian
    % spots and no background noise, so then poissrnd(A_spots_3)
    % generates Poissonian noise for each pixel value (for each pixel, a
    % random value is chosen from a Poissonian distribution of mean equal
    % to the intensity of the pixel value and stdev equal to sqrt of that
    % intensity.

    % Background noise (assume Gaussian (normal) distribution), add for each
    % pixel. The standard deviation of the background noise is bgnd_noise_std:
    bgnd_noise_matrix = offset_level + normrnd(0,bgnd_noise_std,[frame_size frame_size]);
    
    % Include noise on images:
    image_with_noise = A_spots_3 + shot_noise_matrix + bgnd_noise_matrix;
     
    % stack the generated frames onto a frame_size-by-frame_size-by-numFrames
    % matrix, which is the image sequence:
    A_image_sequence(:,:,i) = image_with_noise;
    
end

%% Save final image sequence:

% Save image sequence as a matrix in a .mat file with the corresponding image label:
save(strcat('image_sequence_',image_label),'A_image_sequence')


% Graphic test, display image sequence as a video at the end:
for i=1:numFrames
imshow(A_image_sequence(:,:,i),[]); 
pause(0.2);
end