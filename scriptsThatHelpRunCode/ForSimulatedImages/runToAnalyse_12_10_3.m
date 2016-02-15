
%% Generate simulated images

% CONSTANT PARAMETERS in function simulateImageSequence3.m:

% % For bright spots:
% gaussian_sigma = 3.5; % Gaussian width in pixels (doesn't need to be an integer). (Default = 3.5).
% bleach_spot = 0; % include photobleaching decay of spot intensity (1) or not (0).
% photobleach_tau_spot = 30; % exponential photobleaching decay time constant, in
% % frames (not seconds!). To simulate photobleaching from frame to frame.
% % For time between frames of 40ms, a time constant tau of 1 second
% % corresponds to 25 frames, tau = 2s corresponds to 50 frames, etc.
% % Default: 25-75 frames.  (Default = 30).
% num_bright_spots = 1; % Number of bright spots inside cell autoflu region. (Default = 3).
% num_molecules =1; 
%
% % For cell region, ie., region with autofluorescence:
% autoflu_radius = 20; % radius of the circle mimicking the cell shape and its autofluorescence n.b. 1 extra pixel is added to total diameter by the fspecial function used. (Default = 50).
% bleach_cell = 0; % include photobleaching decay of cell background autoflu intensity (1) or not (0).
% photobleach_tau_autofluo = 25; % rate of bleaching of cell autofluorescence. (Default = 25).
% cell_autofluo_intensity = 0; % mean intensity of  autofluorescence in the simulated cell before photobleaching (scale 0-1). Referred to the amplitude of the bright spot gaussians, which is 1. (Default = 0.3).
% 
% % Noise parameters:
% shot_noise_factor = 0; % For noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity. (Default = 0.5).
% bgnd_noise_std = 0.7 (eg); % For background noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude to noise standard deviation). (Default = 0.3). 
% % Referred to the amplitude of the bright spot gaussians, which is 1.

% PARAMETER TO VARY IS "bgnd_noise_std" in function
% simulateImageSequence3.m.
% In order of images from '01' to '10', "bgnd_noise_std" is: 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.25, 0.35, 0.6, 0.7.  

% A = simulateImageSequence3('01',100,200);
% Bright spot positions x,y(pix): 25    23
% 
% A = simulateImageSequence3('02',100,200);
% Bright spot positions x,y(pix): 25    23
% 
% A = simulateImageSequence3('03',100,200);
% Bright spot positions x,y(pix): 25    27
% 
% A = simulateImageSequence3('04',100,200);
% Bright spot positions x,y(pix):  24    27
% 
% A = simulateImageSequence3('05',100,200);
% Bright spot positions x,y(pix):  23    23
% 
% A = simulateImageSequence3('06',100,200);
% Bright spot positions x,y(pix):  27    23
% 
% A = simulateImageSequence3('07',100,200);
% Bright spot positions x,y(pix):  23    27
% 
% A = simulateImageSequence3('08',100,200);
% Bright spot positions x,y(pix):  24    24
% 
% A = simulateImageSequence3('09',100,200);
% Bright spot positions x,y(pix):  27    27
% 
% A = simulateImageSequence3('10',100,200);
% Bright spot positions x,y(pix):  23    26



%% ANALYSE IMAGE SEQUENCES TO LOOK FOR TRACKS:

% Params for function FindTrajects.m, ~same as for old OXPHOS data:
% r	5
% max_num_candidates	2000
% subarray_halfwidth	10
% inner_circle_radius	8
% gauss_mask_sigma	2
% sigmaFit_min	-8
% sigmaFit_max	8
% SNR_min	2
% rsq_min	0.2
% deflate	1
% d_coincid_cand	1
% d_coincid_found	1
% d_01_max	5
% Iratio_01_min	0.5
% Iratio_01_max	3
% SigmaRatio_01_min	0.5
% SigmaRatio_01_max	2
% d_02_max	5
% Iratio_02_min	0.5
% Iratio_02_max	3
% SigmaRatio_02_min	0.5
% SigmaRatio_02_max	2
% rej	10000
% exclude_region	1

% Params for function linkTrajSegments.m, ~same as for old OXPHOS data:
% rej	10000
% d_01_max	5
% Iratio_01_min	0.2
% Iratio_01_max	4
% SigmaRatio_01_min	0.5
% SigmaRatio_01_max	2
% Frames_away_max	2

% Parameters for function findSpotCentre1frame.m:
% guess_sigma_fit = 3;

% tic;
% s01 = FindTrajects('01',1,200);
% % time elapsed in seconds:
% t_elapsed_01 = toc
% linkTrajSegments('01',1,200,s01,'varyingBgndNoiseStd');
% 
% tic;
% s02 = FindTrajects('02',1,200);
% t_elapsed_02 = toc
% linkTrajSegments('02',1,200,s02,'varyingBgndNoiseStd');
% 
% tic;
% s03 = FindTrajects('03',1,200);
% t_elapsed_03 = toc
% linkTrajSegments('03',1,200,s03,'varyingBgndNoiseStd');
% 
% tic;
% s04 = FindTrajects('04',1,200);
% t_elapsed_04 = toc
% linkTrajSegments('04',1,200,s04,'varyingBgndNoiseStd');


% From here on, change Parameters in FindTrajects.m to:
% SNR_min = 1.4; 
% rsq_min = 0.1;

% tic;
% s05 = FindTrajects('05',1,200);
% t_elapsed_05 = toc
% linkTrajSegments('05',1,200,s05,'varyingBgndNoiseStd');
% 
% tic;
% s06 = FindTrajects('06',1,200);
% t_elapsed_06 = toc
% linkTrajSegments('06',1,200,s06,'varyingBgndNoiseStd');

% tic;
% s07 = FindTrajects('07',1,200);
% t_elapsed_07 = toc
% linkTrajSegments('07',1,200,s07,'varyingBgndNoiseStd');

% tic;
% s08 = FindTrajects('08',1,200);
% t_elapsed_08 = toc
% linkTrajSegments('08',1,200,s08,'varyingBgndNoiseStd');

% tic;
% s09 = FindTrajects('09',1,200);
% t_elapsed_09 = toc
% linkTrajSegments('09',1,200,s09,'varyingBgndNoiseStd');
% 
% tic;
% s10 = FindTrajects('10',1,200);
% t_elapsed_10 = toc
% linkTrajSegments('10',1,200,s10,'varyingBgndNoiseStd');


% --------------------

%% Now adding shot noise to simulated images and adding possibility to
%% modify number of molecules (per bright spot) in image sequence:

% % For bright spots:
% gaussian_sigma = 3.5; % Gaussian width in pixels (doesn't need to be an integer). (Default = 3.5).
% bleach_spot = 0; % include photobleaching decay of spot intensity (1) or not (0).
% photobleach_tau_spot = 30; % exponential photobleaching decay time constant, in
% % frames (not seconds!). To simulate photobleaching from frame to frame.
% % For time between frames of 40ms, a time constant tau of 1 second
% % corresponds to 25 frames, tau = 2s corresponds to 50 frames, etc.
% % Default: 25-75 frames.  (Default = 30).
% num_bright_spots = 1; % Number of bright spots inside cell autoflu region. (Default = 3).
% 
% % For cell region, ie., region with autofluorescence:
% autoflu_radius = 20; % radius of the circle mimicking the cell shape and its autofluorescence n.b. 1 extra pixel is added to total diameter by the fspecial function used. (Default = 50).
% bleach_cell = 0; % include photobleaching decay of cell background autoflu intensity (1) or not (0).
% photobleach_tau_autofluo = 25; % rate of bleaching of cell autofluorescence. (Default = 25).
% cell_autofluo_intensity = 0; % mean intensity of  autofluorescence in the simulated cell before photobleaching (scale 0-1). Referred to the amplitude of the bright spot gaussians, which is 1. (Default = 0.3).
% 
% % Noise parameters:
% shot_noise_factor = 1; % For noise which depends on intensity, for shot-noise: proportional to square root of number of photons/intensity. (Default = 1).
% bgnd_noise_std = 0.5; % For background noise. It is approx. the inverse of the Signal to noise ratio (ratio of signal (gaussian) amplitude to noise standard deviation). (Default = 0.5). 
% % Referred to the amplitude of a single-molecule bright spot Gaussian,
% % which is 1. (obtain from in vitro data)

% Use fixed position of the bright spot for all images, for simplicity:
%     Xpos_bright_spots(k,1) = 25;
%     Ypos_bright_spots(k,1) = 25;

% So use bgnd_noise_std = 0.5, which corresponds to an SNR of around 2,
% which is the approx. one for single-molecule GFP and mCherry from in
% vitro data; assume a shot_noise_factor of 1, and vary the input number of
% molecules, num_molecules, logarithmically (log10) from 1 to 100, in function:
% simulateImageSequence3(image_label,frame_size,numFrames,num_molecules)

% So for noise parameters:
% shot_noise_factor = 1; 
% bgnd_noise_std = 0.5;

% A = simulateImageSequence3('11',100,200,1); 
% 
% A = simulateImageSequence3('12',100,200,2);
% 
% A = simulateImageSequence3('13',100,200,5);
% 
% A = simulateImageSequence3('14',100,200,10);
% 
% A = simulateImageSequence3('15',100,200,20);
% 
% A = simulateImageSequence3('16',100,200,50);
% 
% A = simulateImageSequence3('17',100,200,100);




% Same but for noise parameters:
% shot_noise_factor = 1; 
% bgnd_noise_std = 0.6;

% A = simulateImageSequence3('21',100,200,1);  
% 
% A = simulateImageSequence3('22',100,200,2);
% 
% A = simulateImageSequence3('23',100,200,5);
% 
% A = simulateImageSequence3('24',100,200,10);
% 
% A = simulateImageSequence3('25',100,200,20);
% 
% A = simulateImageSequence3('26',100,200,50);
% 
% A = simulateImageSequence3('27',100,200,100);




% Same but for noise parameters:
% shot_noise_factor = 1; 
% bgnd_noise_std = 0.4;

% A = simulateImageSequence3('31',100,200,1);
% 
% A = simulateImageSequence3('32',100,200,2);
% 
% A = simulateImageSequence3('33',100,200,5);
% 
% A = simulateImageSequence3('34',100,200,10);
% 
% A = simulateImageSequence3('35',100,200,20);
% 
% A = simulateImageSequence3('36',100,200,50);
% 
% A = simulateImageSequence3('37',100,200,100);



% Same but for noise parameters:
% shot_noise_factor = 0; 
% bgnd_noise_std = 0.5;

% A = simulateImageSequence3('41',100,200,1);
% 
% A = simulateImageSequence3('42',100,200,2);
% 
% A = simulateImageSequence3('43',100,200,5);
% 
% A = simulateImageSequence3('44',100,200,10);
% 
% A = simulateImageSequence3('45',100,200,20);
% 
% A = simulateImageSequence3('46',100,200,50);
% 
% A = simulateImageSequence3('47',100,200,100);


% Analysis, tracking of previous images 11 to 47:

% s11 = FindTrajects('11',1,200);
% linkTrajSegments('11',1,200,s11,'varyNumMolecs');
% 
% s12 = FindTrajects('12',1,200);
% linkTrajSegments('12',1,200,s12,'varyNumMolecs');
% 
% s13 = FindTrajects('13',1,200);
% linkTrajSegments('13',1,200,s13,'varyNumMolecs');
% 
% s14 = FindTrajects('14',1,200);
% linkTrajSegments('14',1,200,s14,'varyNumMolecs');
% 
% s15 = FindTrajects('15',1,200);
% linkTrajSegments('15',1,200,s15,'varyNumMolecs');
% 
% s16 = FindTrajects('16',1,200);
% linkTrajSegments('16',1,200,s16,'varyNumMolecs');
% 
% s17 = FindTrajects('17',1,200);
% linkTrajSegments('17',1,200,s17,'varyNumMolecs');
% 
% 
% 
% s21 = FindTrajects('21',1,200);
% linkTrajSegments('21',1,200,s21,'varyNumMolecs');
% 
% s22 = FindTrajects('22',1,200);
% linkTrajSegments('22',1,200,s22,'varyNumMolecs');
% 
% s23 = FindTrajects('23',1,200);
% linkTrajSegments('23',1,200,s23,'varyNumMolecs');
% 
% s24 = FindTrajects('24',1,200);
% linkTrajSegments('24',1,200,s24,'varyNumMolecs');
% 
% s25 = FindTrajects('25',1,200);
% linkTrajSegments('25',1,200,s25,'varyNumMolecs');
% 
% s26 = FindTrajects('26',1,200);
% linkTrajSegments('26',1,200,s26,'varyNumMolecs');
% 
% s27 = FindTrajects('27',1,200);
% linkTrajSegments('27',1,200,s27,'varyNumMolecs');
% 
% 
% 
% s31 = FindTrajects('31',1,200);
% linkTrajSegments('31',1,200,s31,'varyNumMolecs');
% 
% s32 = FindTrajects('32',1,200);
% linkTrajSegments('32',1,200,s32,'varyNumMolecs');
% 
% s33 = FindTrajects('33',1,200);
% linkTrajSegments('33',1,200,s33,'varyNumMolecs');
% 
% s34 = FindTrajects('34',1,200);
% linkTrajSegments('34',1,200,s34,'varyNumMolecs');
% 
% s35 = FindTrajects('35',1,200);
% linkTrajSegments('35',1,200,s35,'varyNumMolecs');
% 
% s36 = FindTrajects('36',1,200);
% linkTrajSegments('36',1,200,s36,'varyNumMolecs');
% 
% s37 = FindTrajects('37',1,200);
% linkTrajSegments('37',1,200,s37,'varyNumMolecs');
% 
% 
% 
% s41 = FindTrajects('41',1,200);
% linkTrajSegments('41',1,200,s41,'varyNumMolecs');
% 
% s42 = FindTrajects('42',1,200);
% linkTrajSegments('42',1,200,s42,'varyNumMolecs');
% 
% s43 = FindTrajects('43',1,200);
% linkTrajSegments('43',1,200,s43,'varyNumMolecs');
% 
% s44 = FindTrajects('44',1,200);
% linkTrajSegments('44',1,200,s44,'varyNumMolecs');
% 
% s45 = FindTrajects('45',1,200);
% linkTrajSegments('45',1,200,s45,'varyNumMolecs');
% 
% s46 = FindTrajects('46',1,200);
% linkTrajSegments('46',1,200,s46,'varyNumMolecs');
% 
% s47 = FindTrajects('47',1,200);
% linkTrajSegments('47',1,200,s47,'varyNumMolecs');


%% Error analysis: compare truth to results from tracking-spotFinding;

% % shot_noise_factor = 1:
% r11 = compareToTruth('11',25,25,3.5,0,0.5,1);
% r12 = compareToTruth('12',25,25,3.5,0,0.5,2);
% r13 = compareToTruth('13',25,25,3.5,0,0.5,5);
% r14 = compareToTruth('14',25,25,3.5,0,0.5,10);
% r15 = compareToTruth('15',25,25,3.5,0,0.5,20);
% r16 = compareToTruth('16',25,25,3.5,0,0.5,50);
% r17 = compareToTruth('17',25,25,3.5,0,0.5,100);
% 
% r21 = compareToTruth('21',25,25,3.5,0,0.6,1);
% r22 = compareToTruth('22',25,25,3.5,0,0.6,2);
% r23 = compareToTruth('23',25,25,3.5,0,0.6,5);
% r24 = compareToTruth('24',25,25,3.5,0,0.6,10);
% r25 = compareToTruth('25',25,25,3.5,0,0.6,20);
% r26 = compareToTruth('26',25,25,3.5,0,0.6,50);
% r27 = compareToTruth('27',25,25,3.5,0,0.6,100);
% 
% r31 = compareToTruth('31',25,25,3.5,0,0.4,1);
% r32 = compareToTruth('32',25,25,3.5,0,0.4,2);
% r33 = compareToTruth('33',25,25,3.5,0,0.4,5);
% r34 = compareToTruth('34',25,25,3.5,0,0.4,10);
% r35 = compareToTruth('35',25,25,3.5,0,0.4,20);
% r36 = compareToTruth('36',25,25,3.5,0,0.4,50);
% r37 = compareToTruth('37',25,25,3.5,0,0.4,100);
% 
% % shot_noise_factor = 0:
% r41 = compareToTruth('41',25,25,3.5,0,0.5,1);
% r42 = compareToTruth('42',25,25,3.5,0,0.5,2);
% r43 = compareToTruth('43',25,25,3.5,0,0.5,5);
% r44 = compareToTruth('44',25,25,3.5,0,0.5,10);
% r45 = compareToTruth('45',25,25,3.5,0,0.5,20);
% r46 = compareToTruth('46',25,25,3.5,0,0.5,50);
% r47 = compareToTruth('47',25,25,3.5,0,0.5,100);
% 
% save 'meanErrorsComparedToTruth' 'r*'


%% Plot results:

load 'meanErrorsComparedToTruth'

numMolecs_axis = [1 2 5 10 20 50 100]'; % column vector.

results_vector_1 = [r11 r12 r13 r14 r15 r16 r17];
results_vector_2 = [r21 r22 r23 r24 r25 r26 r27];
results_vector_3 = [r31 r32 r33 r34 r35 r36 r37];
results_vector_4 = [r41 r42 r43 r44 r45 r46 r47];

field_names = fields(r11); % cell structure
% 1    'localis_errorX_pix_abs'
% 2   'localis_errorY_pix_abs'
% 3    'localis_errorX_percent_abs'
% 4    'localis_errorY_percent_abs'
% 5    'Sigma_error_percent_abs'
% 6    'IspTot_error_percent_abs'
% 7    'NumMolecs_error_percent_abs'
% 8     'OffsetBgNoise_error_abs'
% 9    'BgNoiseStd_error_percent_abs'
% 10    'SNR_error_percent_abs'
% 11   'localis_errorX_pix'
% 12  'localis_errorY_pix'
% 13  'localis_errorX_percent'
% 14  'localis_errorY_percent'
% 15  'Sigma_error_percent'
% 16  'IspTot_error_percent'
% 17  'NumMolecs_error_percent'
% 18    'OffsetBgNoise_error'
% 19  'BgNoiseStd_error_percent'
% 20  'SNR_error_percent'

% The output of function "getErrorVectors" is a matrix in which each column
% is a variable (each field within the result structures r11, r12 etc...), and
% each row corresponds to a given number of molecules:
yVectorsToPlot_1 = getErrorVectors(results_vector_1,numMolecs_axis);
yVectorsToPlot_2 = getErrorVectors(results_vector_2,numMolecs_axis);
yVectorsToPlot_3 = getErrorVectors(results_vector_3,numMolecs_axis);
yVectorsToPlot_4 = getErrorVectors(results_vector_4,numMolecs_axis);

% % Plot localisation errors in pixels:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,1),...
%     numMolecs_axis,yVectorsToPlot_3(:,2),...
%     numMolecs_axis,yVectorsToPlot_1(:,1),...
%     numMolecs_axis,yVectorsToPlot_1(:,2),...
%     numMolecs_axis,yVectorsToPlot_2(:,1),...
%     numMolecs_axis,yVectorsToPlot_2(:,2),...
%     numMolecs_axis,yVectorsToPlot_4(:,1),...
%     numMolecs_axis,yVectorsToPlot_4(:,2));
% xlabel('number of molecules');
% ylabel('localisation error (pixels)');
% legend('X, bgNoiseStd=0.4, with shot noise','Y, bgNoiseStd=0.4, with shot noise',...
%     'X, bgNoiseStd=0.5, with shot noise','Y, bgNoiseStd=0.5, with shot noise',...
%     'X, bgNoiseStd=0.6, with shot noise','Y, bgNoiseStd=0.6, with shot noise',...
%     'X, bgNoiseStd=0.5, no shot noise','Y, bgNoiseStd=0.5, no shot noise')
% 
% % Plot localisation errors in percentage:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,3),...
%     numMolecs_axis,yVectorsToPlot_3(:,4),...
%     numMolecs_axis,yVectorsToPlot_1(:,3),...
%     numMolecs_axis,yVectorsToPlot_1(:,4),...
%     numMolecs_axis,yVectorsToPlot_2(:,3),...
%     numMolecs_axis,yVectorsToPlot_2(:,4),...
%     numMolecs_axis,yVectorsToPlot_4(:,3),...
%     numMolecs_axis,yVectorsToPlot_4(:,4));
% xlabel('number of molecules');
% ylabel('localisation error (%)');
% legend('X, bgNoiseStd=0.4, with shot noise','Y, bgNoiseStd=0.4, with shot noise',...
%     'X, bgNoiseStd=0.5, with shot noise','Y, bgNoiseStd=0.5, with shot noise',...
%     'X, bgNoiseStd=0.6, with shot noise','Y, bgNoiseStd=0.6, with shot noise',...
%     'X, bgNoiseStd=0.5, no shot noise','Y, bgNoiseStd=0.5, no shotnoise')
%     
% 
% % Plot localisation errors in pixels:
% % sqrt(errorX^2+errorY^2):
% figure
% loglog(numMolecs_axis,sqrt(yVectorsToPlot_3(:,1).^2+yVectorsToPlot_3(:,2).^2),...
%     numMolecs_axis,sqrt(yVectorsToPlot_1(:,1).^2+yVectorsToPlot_1(:,2).^2),...
%     numMolecs_axis,sqrt(yVectorsToPlot_2(:,1).^2+yVectorsToPlot_2(:,2).^2),...
%     numMolecs_axis,sqrt(yVectorsToPlot_4(:,1).^2+yVectorsToPlot_4(:,2).^2));
% xlabel('number of molecules');
% ylabel('Localisation error (pixels)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise');
% 
% % Plot localisation errors in nanometres:
% % sqrt(errorX^2+errorY^2):
% pixelToNM = 35.55; % conversion from pixels to nm.
% figure
% loglog(numMolecs_axis,pixelToNM*sqrt(yVectorsToPlot_3(:,1).^2+yVectorsToPlot_3(:,2).^2),...
%     numMolecs_axis,pixelToNM*sqrt(yVectorsToPlot_1(:,1).^2+yVectorsToPlot_1(:,2).^2),...
%     numMolecs_axis,pixelToNM*sqrt(yVectorsToPlot_2(:,1).^2+yVectorsToPlot_2(:,2).^2),...
%     numMolecs_axis,pixelToNM*sqrt(yVectorsToPlot_4(:,1).^2+yVectorsToPlot_4(:,2).^2));
% xlabel('number of molecules');
% ylabel('Localisation error (nm)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise');
% 
% % Plot localisation errors in pixels:
% % Only for x:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,1),...
%     numMolecs_axis,yVectorsToPlot_1(:,1),...
%     numMolecs_axis,yVectorsToPlot_2(:,1),...
%     numMolecs_axis,yVectorsToPlot_4(:,1));
% xlabel('number of molecules');
% ylabel('X localisation error (pixels)');
% legend('X, bgNoiseStd=0.4, with shot noise',...
%     'X, bgNoiseStd=0.5, with shot noise',...
%     'X, bgNoiseStd=0.6, with shot noise',...
%     'X, bgNoiseStd=0.5, no shot noise');
% 
% % Only for y:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,2),...
%     numMolecs_axis,yVectorsToPlot_1(:,2),...
%     numMolecs_axis,yVectorsToPlot_2(:,2),...
%     numMolecs_axis,yVectorsToPlot_4(:,2));
% xlabel('number of molecules');
% ylabel('Y localisation error (pixels)');
% legend('Y, bgNoiseStd=0.4, with shot noise',...
%     'Y, bgNoiseStd=0.5, with shot noise',...
%     'Y, bgNoiseStd=0.6, with shot noise',...
%     'Y, bgNoiseStd=0.5, no shot noise');
% xlim([0 10]) % zoom in for x axis.
% 
% 
% % Plot percentage error in width:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,5),...
%     numMolecs_axis,yVectorsToPlot_1(:,5),...
%     numMolecs_axis,yVectorsToPlot_2(:,5),...
%     numMolecs_axis,yVectorsToPlot_4(:,5));
% xlabel('number of molecules');
% ylabel('Width error (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in number of molecues (or I0):
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,7),...
%     numMolecs_axis,yVectorsToPlot_1(:,7),...
%     numMolecs_axis,yVectorsToPlot_2(:,7),...
%     numMolecs_axis,yVectorsToPlot_4(:,7));
% xlabel('number of molecules');
% ylabel('NumMolecs error (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot absolute error in Offset Background level after bgnd subtraction (compared to 0):
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,8),...
%     numMolecs_axis,yVectorsToPlot_1(:,8),...
%     numMolecs_axis,yVectorsToPlot_2(:,8),...
%     numMolecs_axis,yVectorsToPlot_4(:,8));
% xlabel('number of molecules');
% ylabel('Offset background level (avg counts per pixel)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in bg_noise_std:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,9),...
%     numMolecs_axis,yVectorsToPlot_1(:,9),...
%     numMolecs_axis,yVectorsToPlot_2(:,9),...
%     numMolecs_axis,yVectorsToPlot_4(:,9));
% xlabel('number of molecules');
% ylabel('bg-noise-stdev error (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in SNR:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,10),...
%     numMolecs_axis,yVectorsToPlot_1(:,10),...
%     numMolecs_axis,yVectorsToPlot_2(:,10),...
%     numMolecs_axis,yVectorsToPlot_4(:,10));
% xlabel('number of molecules');
% ylabel('SNR error (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in IspTot:
% figure
% loglog(numMolecs_axis,yVectorsToPlot_3(:,6),...
%     numMolecs_axis,yVectorsToPlot_1(:,6),...
%     numMolecs_axis,yVectorsToPlot_2(:,6),...
%     numMolecs_axis,yVectorsToPlot_4(:,6));
% xlabel('number of molecules');
% ylabel('IspTot error (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')




%% Now the non-absolute value ones, in case things show a systematic
%% error in one direction:

% % Plot localisation errors in percentage:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,13),...
%     numMolecs_axis,yVectorsToPlot_3(:,14),...
%     numMolecs_axis,yVectorsToPlot_1(:,13),...
%     numMolecs_axis,yVectorsToPlot_1(:,14),...
%     numMolecs_axis,yVectorsToPlot_2(:,13),...
%     numMolecs_axis,yVectorsToPlot_2(:,14),...
%     numMolecs_axis,yVectorsToPlot_4(:,13),...
%     numMolecs_axis,yVectorsToPlot_4(:,14));
% xlabel('number of molecules');
% ylabel('localisation error with sign (%)');
% legend('X, bgNoiseStd=0.4, with shot noise','Y, bgNoiseStd=0.4, with shot noise',...
%     'X, bgNoiseStd=0.5, with shot noise','Y, bgNoiseStd=0.5, with shot noise',...
%     'X, bgNoiseStd=0.6, with shot noise','Y, bgNoiseStd=0.6, with shot noise',...
%     'X, bgNoiseStd=0.5, no shot noise','Y, bgNoiseStd=0.5, no shotnoise')
% 
% 
% % Plot percentage error in width:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,15),...
%     numMolecs_axis,yVectorsToPlot_1(:,15),...
%     numMolecs_axis,yVectorsToPlot_2(:,15),...
%     numMolecs_axis,yVectorsToPlot_4(:,15));
% xlabel('number of molecules');
% ylabel('Width error with sign (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in number of molecues (or I0):
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,17),...
%     numMolecs_axis,yVectorsToPlot_1(:,17),...
%     numMolecs_axis,yVectorsToPlot_2(:,17),...
%     numMolecs_axis,yVectorsToPlot_4(:,17));
% xlabel('number of molecules');
% ylabel('NumMolecs error with sign (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot error in offset background noise as difference with real 0 offset:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,18),...
%     numMolecs_axis,yVectorsToPlot_1(:,18),...
%     numMolecs_axis,yVectorsToPlot_2(:,18),...
%     numMolecs_axis,yVectorsToPlot_4(:,18));
% xlabel('number of molecules');
% ylabel('offset noise compared to 0 real offset');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in bg_noise_std:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,19),...
%     numMolecs_axis,yVectorsToPlot_1(:,19),...
%     numMolecs_axis,yVectorsToPlot_2(:,19),...
%     numMolecs_axis,yVectorsToPlot_4(:,19));
% xlabel('number of molecules');
% ylabel('bg-noise-stdev error with sign (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in SNR:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,20),...
%     numMolecs_axis,yVectorsToPlot_1(:,20),...
%     numMolecs_axis,yVectorsToPlot_2(:,20),...
%     numMolecs_axis,yVectorsToPlot_4(:,20));
% xlabel('number of molecules');
% ylabel('SNR error with sign (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% 
% 
% % Plot percentage error in IspTot:
% figure
% plot(numMolecs_axis,yVectorsToPlot_3(:,16),...
%     numMolecs_axis,yVectorsToPlot_1(:,16),...
%     numMolecs_axis,yVectorsToPlot_2(:,16),...
%     numMolecs_axis,yVectorsToPlot_4(:,16));
% xlabel('number of molecules');
% ylabel('IspTot error with sign (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')


%% Percentage of spots detected:
%(length of excel file/2), out of 200 frames.

% In increasing order of number of molecules (1,2,5,10,20,50,100), and as a percentage:
fraction_spotsDetected_1 = [42,98,100,100,100,100,100];
fraction_spotsDetected_2 = [18,100,100,100,100,100,100];
fraction_spotsDetected_3 = [58,100,100,100,100,100,100];
fraction_spotsDetected_4 = [99,100,100,100,100,100,100];

% figure
% loglog(numMolecs_axis,fraction_spotsDetected_3,...
%     numMolecs_axis,fraction_spotsDetected_1,...
%     numMolecs_axis,fraction_spotsDetected_2,...
%     numMolecs_axis,fraction_spotsDetected_4);
% xlabel('number of molecules');
% ylabel('Fraction of spots detected (%)');
% legend('bgNoiseStd=0.4, with shot noise',...
%     'bgNoiseStd=0.5, with shot noise',...
%     'bgNoiseStd=0.6, with shot noise',...
%     'bgNoiseStd=0.5, no shot noise')
% ylim([10 120])
