function A_image_sequence = simulateImageSequence3(image_label,frame_size,numFrames,num_molecules)
%
% Started by Isabel Llorente-Garcia. April 2012.
% Modified and expanded by George Averill in Sept. 2012.
% Then modified again by I. Ll.-G. October 2012.
%
% Similar to simulateImageSequence.m but adding a background of cell
% autofluorescence. There is the lowest background, then a circular region
% of higher cell background autofluorescence and then, on the top channel
% only, smaller Gaussian bright spots within the autofluorescent region.
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
% - num_molecules % number of single-molecules (Gaussian intensity spots with amplitude I0 equal to 1) within one bright spot.

%
% OUTPUTS:
% A .mat file is saved with a MxMxJ matrix, where MxM is the size of each
% frame (M=frame_size), and J=numFrames is the number of frames in the image sequence.
%
% EXAMPLE of how to call this function: 
% A = simulateImageSequence3('try1',100,200,5);
%
% Note: gaussian bright spots have an amplitude of 1, and cell-autoflu has
% an additional amplitude of cell_autofluo_intensity. And noise added on
% top of that. 


%% PARAMETERS:

% For bright spots:
gaussian_sigma = 3.5; % Gaussian width in pixels (doesn't need to be an integer). (Default = 3.5).
bleach_spot = 0; % include photobleaching decay of spot intensity (1) or not (0).
photobleach_tau_spot = 30; % exponential photobleaching decay time constant, in
% frames (not seconds!). To simulate photobleaching from frame to frame.
% For time between frames of 40ms, a time constant tau of 1 second
% corresponds to 25 frames, tau = 2s corresponds to 50 frames, etc.
% Default: 25-75 frames.  (Default = 30).
num_bright_spots = 1; % Number of bright spots inside cell autoflu region. (Default = 3).

% For cell region, ie., region with autofluorescence:
autoflu_radius = 20; % radius of the circle mimicking the cell shape and its autofluorescence n.b. 1 extra pixel is added to total diameter by the fspecial function used. (Default = 50).
bleach_cell = 0; % include photobleaching decay of cell background autoflu intensity (1) or not (0).
photobleach_tau_autofluo = 25; % rate of bleaching of cell autofluorescence. (Default = 25).
cell_autofluo_intensity = 0; % mean intensity of  autofluorescence in the simulated cell before photobleaching (scale 0-1). Referred to the amplitude of the bright spot gaussians, which is 1. (Default = 0.3).

% Noise parameters:
shot_noise_factor = 1; % For noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity. (Default = 1).
bgnd_noise_std = 0.5; % For background noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude to noise standard deviation). (Default = 0.5). 
% Referred to the amplitude of a single-molecule bright spot Gaussian,
% which is 1. (obtain from in vitro data)


%% Error control:

if round(frame_size/4) < autoflu_radius
    disp ('----  !!!!!!  ----')
    disp('ERROR!! - Check parameters: autoflu_radius needs to be <= frame_size/4 !!!')
    disp('Exiting function.')
    disp ('----  !!!!!!  ----')
    return
end

%% Generate simulated image sequence:

% Initialise final image-sequence matrix to store things in:
A_image_sequence = zeros(frame_size,frame_size,numFrames);

% SIMULATED CIRCULAR CELL:
% Initialise top-left quarter of image with simulated cell, for first frame:
A_top_left_cell_1 = zeros(frame_size/2,frame_size/2);
% Create uniform circle to simulate the autofluorescence from the cell region (simulate a circular cell):
circle = fspecial('disk', autoflu_radius); % Creates uniform intensity circle of radius 'autoflu_radius' on a black square background
% Position for simulated cell centre: pixels to locate simulated circular cell region:
x_cell_centre = round(frame_size/4); % centre of top-left quarter of image.
y_cell_centre = round(frame_size/4); % centre of top-left quarter of image.
A_top_left_cell_1(y_cell_centre,x_cell_centre) = 1; % place a 1 at the position of the cell centre in A_top_left_cell_1.
% Filter with "circle" to get the simulated circular cell at the desired position:
A_top_left_cell_2 = imfilter(A_top_left_cell_1,circle);
A_top_left_cell_3 = mat2gray(A_top_left_cell_2)*cell_autofluo_intensity; % scales between 0 and 1 (to greyscale image) and then scales for the intensity parameter.

% SIMULATED BRIGHT FLUORESCENT SPOTS:
% Generate a small matrix with one Gaussian bright spot in it (to act as filter later on, to simulate flu bright spots):
gauss_bright_spot = fspecial('gaussian',min(round(10*gaussian_sigma),frame_size),gaussian_sigma); % square matrix of size 10*gaussian_sigma, with a normalised Gaussian spot of sigma "gaussian_sigma" pixels. 
% Generate random positions for bright spots inside top-left quarter of image:
Xpos_bright_spots = zeros(num_bright_spots,1); % initialise column vector.
Ypos_bright_spots = zeros(num_bright_spots,1); % initialise column vector.
% Initialise top-left quarter of image with bright spots, for first frame:
A_top_left_spots_1 = zeros(frame_size/2,frame_size/2);
for k=1:num_bright_spots
    % Make sure the centre of the bright spot lies well inside the simulated
    % cell region:   
%     Xpos_bright_spots(k,1) = randi([x_cell_centre-autoflu_radius+round(5*gaussian_sigma),x_cell_centre+autoflu_radius-round(5*gaussian_sigma)]);
%     Ypos_bright_spots(k,1) = randi([y_cell_centre-autoflu_radius+round(5*gaussian_sigma),y_cell_centre+autoflu_radius-round(5*gaussian_sigma)]);
    % Or use fixed positions of choice  (use for a single bright spot only: num_bright_spots = 1!)(and comment out previous two lines):
    Xpos_bright_spots(k,1) = 25;
    Ypos_bright_spots(k,1) = 25;
%     % Make sure the centre of the bright spot lies inside the simulated cell region:
%     while A_top_left_cell_3(Ypos_bright_spots(k,1),Xpos_bright_spots(k,1)) < cell_autofluo_intensity
%             Xpos_bright_spots(k,1) = randi(frame_size/2);
%             Ypos_bright_spots(k,1) = randi(frame_size/2);
%     end
    A_top_left_spots_1(Ypos_bright_spots(k,1),Xpos_bright_spots(k,1)) = 1; % place a 1 at the positions of bright spots in A_top_left_quarter_1b.
end

disp('Bright spot positions x,y(pix):')
disp([Xpos_bright_spots Ypos_bright_spots])

% Place gaussian bright flu spots at positions given by above places, use filter:
A_top_left_spots_2 = imfilter(A_top_left_spots_1,gauss_bright_spot);
% Scale between 0 and 1 (to grayscale image), and multiply by the number of molecules per bright spot (each spot has a maximum Gaussian amplitude of 1):
A_top_left_spots_3 = num_molecules*mat2gray(A_top_left_spots_2); 
% Asume I0, the amplitude of a Gaussian intensity spot from a single
% molecule is 1. Then assume width stays the same when increasing the
% number of molecules and that I0 scales with the number of molecules,
% num_molecules, in one bright spot.


% LOOP THROUGH FRAMES:
for i=1:numFrames
    
    % Account for exponential photobleaching decay of bright spots:
    if bleach_spot == 1
        bleach_factor_spots = exp(-(i-1)/photobleach_tau_spot);
    else
        bleach_factor_spots = 1;
    end
    A_top_left_spots_4 = bleach_factor_spots*A_top_left_spots_3;
    
    % Account for exponential photobleaching decay of cell
    % autofluorescence:
    if bleach_cell == 1
        bleach_factor_cell_autoflu = exp(-(i-1)/photobleach_tau_autofluo);
    else
        bleach_factor_cell_autoflu = 1;
    end
    A_top_left_cell_4 = bleach_factor_cell_autoflu*A_top_left_cell_3;
    
    % Overlay noiseless cell autoflu region and bright spots (top left quadrant of image): 
    A_top_left_quadrant = A_top_left_cell_4 + A_top_left_spots_4;
    % Note that this is a square matrix of size frame_size/2, for which the
    % max value is 1+cell_autofluo_intensity.
    
    padded_image = padarray(A_top_left_quadrant,[frame_size/2 frame_size/2], 'post'); % adds the other image quadrants to the right and bottom.
        
    % Add noise to images:
    % Shot-noise: proportional to square root of total number of photons/intensity:
    % The uncertainty in the no. of counts (standard deviation of the
    % noise) is proportional to the square root of the number of counts at each position:
    shot_noise_matrix = normrnd(0,shot_noise_factor*sqrt(padded_image));
    % Note, this might be more easily (and correctly) done using the
    % function poissrnd. So if padded_image is the image with my Gaussian
    % spots and no background noise, I can then do poissrnd(padded_image)
    % to generate Poissonian noise for each pixel value (for each pixel, a
    % random value is chosen from a Poissonian distribution of mean equal
    % to the intensity of the pixel value and stdev equal to sqrt of that
    % intensity.


    % Background noise:
    % The standard deviation of the background noise is bgnd_noise_std:
    bgnd_noise_matrix = normrnd(0,bgnd_noise_std,[frame_size frame_size]);
    
    % Include noise on images:
    padded_image_with_noise = padded_image + shot_noise_matrix + bgnd_noise_matrix;
     
    % stack the generated frames onto a frame_size-by-frame_size-by-numFrames
    % matrix, which is the image sequence:
    A_image_sequence(:,:,i) = padded_image_with_noise;
    
end

%% Save final image sequence:

% Save image sequence as a matrix in a .mat file with the corresponding image label:
save(strcat('image_sequence_',image_label),'A_image_sequence')


% Graphic test, display image sequence as a video at the end:
for i=1:numFrames
imshow(A_image_sequence(:,:,i),[]); 
pause(0.2);
end


