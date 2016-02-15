% 
% Run this file to analyse the in vitro data (mCherry fixed to slide with antibodies).
% The ROI is entered by hand and saved onto the results as parameters.
% A different folder is created for the results for each ROI.
%
% 


%% Frame averages:
% --------------------------------------------
% % frameAverage(image_label,start_frame,end_frame,display_result,save_png)

% frameAverage('1',1,100,1,1);
% frameAverage('2',1,100,1,1);
% frameAverage('3',1,100,1,1);
% frameAverage('4',1,100,1,1);
% frameAverage('5good',1,200,1,1);
% frameAverage('6',1,200,1,1);
% frameAverage('7bright',1,200,1,1);
% frameAverage('8',1,200,1,1);
% frameAverage('9',1,200,1,1);
% frameAverage('10good',1,200,1,1);
% frameAverage('11good',1,200,1,1);
% --------------



%% Beam profile: intensity variation accross the image plane:
% -----------------------------------------------------------
% See samplingRowsAndColumnsInFrameAvg.m 
% See saved images in folder "profileBeamIntensity".
% --------------


%% Analyse in vitro bright spots:

% Reminder: 
% result_final = analyseFixedBrightSpots(image_label,start_frame,end_frame,folder_label)

% ---------------------------------------------------------
% PARAMS, eg.:
%                     image_label: '1'
%                    folder_label: 'mCherry'
%                     start_frame: 20
%                       end_frame: 100
%              max_num_candidates: 10000
%              subarray_halfwidth: 8
%             inner_circle_radius: 5
%                gauss_mask_sigma: 2
%                 guess_sigma_Fit: 3
%                    sigmaFit_min: -5
%                    sigmaFit_max: 5
%                         SNR_min: 1.6000
%                         rsq_min: 0.2000
%                         deflate: 1
%                 square_halfsize: 10
%                         Wfilter: 3
%                         Rfilter: 1
%                           nbins: 100
%                  x_limit_spectr: 5000
%          subtract_bgnd_fit_hist: 0
%                    use_filtered: 0
%                 doPairWiseDiffs: 1
%                      image_path: '1.sif'
%                       numFrames: 100
%                     frame_Ysize: 512
%                     frame_Xsize: 512
%                       xleft_ROI: 225
%                      xright_ROI: 315
%                        ytop_ROI: 260
%                     ybottom_ROI: 350
%                     bgnd_spot_x: 20
%                     bgnd_spot_y: 50
%                           Ibgnd: [100x1 double]
%                  Ibgnd_filtered: [100x1 double]
%       Ibgnd_filtered_normalised: [100x1 double]
%                binCentres_Ibgnd: [100x1 double]
%                freqCounts_Ibgnd: [100x1 double]
%           power_spectrum_X_bgnd: [65x1 double]
%           power_spectrum_Y_bgnd: [65x1 double]
%     power_spectrum_peaks_X_bgnd: [17x1 double]
%     power_spectrum_peaks_Y_bgnd: [17x1 double]
%                  Ibgnd_rmsNoise: 6.9762e+003
%                 Ibgnd_ampli_fit: 48.8303
% ---------------------------------------------------------

% Go through images:
% ---------------------------------------------------------
% Set  ROI = [1 512 1 512];
% and y-axis scaling for zoom out plots to limitsYzoom = [-2000 10000];.

% r1 = analyseFixedBrightSpots('1',20,100,'');
% 
% r2 = analyseFixedBrightSpots('2',20,100,'');
% 
% r3 = analyseFixedBrightSpots('3',20,100,'');
% 
% r4 = analyseFixedBrightSpots('4',20,100,'');
% 
% r5 = analyseFixedBrightSpots('5good',20,100,'');
% 
% r6 = analyseFixedBrightSpots('6',30,100,'');
% 
% r7 = analyseFixedBrightSpots('7bright',20,100,'');
% 
% r8 = analyseFixedBrightSpots('8',20,100,'');
% 
% r9 = analyseFixedBrightSpots('9',20,100,'');
% 
% r10 = analyseFixedBrightSpots('10good',20,100,'');
% 
% r11 = analyseFixedBrightSpots('11good',20,100,'');


% 
% Save results:
% ---------------------------------------------------------
% save 'results' 'r*' % save all results in a .mat file in current folder.



%% Put together intensities of all selected spots and make histogram of them, with threshold:
% To try and distinguish different intensity levels.

load results.mat;

% NOTE: to use this first comment out the second part ("Method 2" section)
% of function "analyseIstepResults.m".
disp(' '); % empty line
disp('For ALL spots:');
allImages = {r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11}; % results from all images.
intensities_vector_allImages = []; % initialise vector to accummulate intensity values from all images;
threshold = -Inf; % only accummulate intensity values above this threshold (use -Inf to use all values).   
for j =1:length(allImages)
    
    r_toAnalyse = allImages{j};
    
    r_Ivectors_0 = [r_toAnalyse{1}.integrated_Ispot_bgndSubtracted]; % matrix in which each column is the intensity of a given spot, each row is frame no.
    r_Ivectors = r_Ivectors_0(21:end,:); % eliminate the first 20 frames (all values are 0 because these frames were not analysed). (otherwise we get a very large peak at zero.)
    intensities_vector = []; % initialise vector to accummulate intensity values;
    for i = 1:size(r_Ivectors,2)
        % go through selected spots in these image:
        intensity_values = r_Ivectors(:,i); % column vector.
        intensities_vector = [intensities_vector; intensity_values(intensity_values > threshold)]; % append values larger than threshold;
    end
    
    intensities_vector_allImages = [intensities_vector_allImages; intensities_vector];
    
end
% % Results are robust to change in the number of histogram bins, as they should:
% analyseIstepResults(intensities_vector_allImages,250,1400,250,1400,1)
% analyseIstepResults(intensities_vector_allImages,500,1400,500,1400,1)
% analyseIstepResults(intensities_vector_allImages,1000,1400,1000,1400,1)

% % Check also all the spots corresponding to a single image sequence (vector intensities_vector) (put a breakpoint inside loop):
% % Calculate and plot histogram of all intensity values and also of their pair-wise differences:
% analyseIstepResults(intensities_vector,100,4000,100,4000,1) % make the x-axis max limit of the Fourier spectra equal to 4000.
% analyseIstepResults(intensities_vector,100,2000,100,2000,1) % make the x-axis max limit of the Fourier spectra equal to 2000 (zoom in).
% analyseIstepResults(intensities_vector,100,1000,100,1000,1) % make the x-axis max limit of the Fourier spectra equal to 1000 (zoom further in).


% Save to existing output excel file: 
% CAREFUL!!! It will overwrite existing file!
xlswrite('mCherry_invitroResults',intensities_vector_allImages,'all_I_values'); % save to sheet 'all_I_values_noThreshold'.
% Calculate and plot histogram of all values:
% x and y values of the histogram plot:
[freqCounts,binCentres] = hist(intensities_vector_allImages,250);
[freqCounts1,binCentres1] = hist(intensities_vector_allImages,500);
[freqCounts2,binCentres2] = hist(intensities_vector_allImages,1000);
% Save to existing output excel file: 
% CAREFUL!!! It will overwrite existing file!
xlswrite('mCherry_invitroResults',[binCentres',freqCounts'],'histogram x-y all_250bars');
xlswrite('mCherry_invitroResults',[binCentres1',freqCounts1'],'histogram x-y all_500bars');
xlswrite('mCherry_invitroResults',[binCentres2',freqCounts2'],'histogram x-y all_1000bars');

% Plot histogram in green:
bar(binCentres2,freqCounts2,'g');
hold on;

% Take only the left hand side of main peak (highlighted in blue on plot) and only up to approx. the FWHM to its right hand side and fit it to a gaussian.
% We assume the distribution is Gaussian. Fit shown as a black curve.
x_lim_for_fit = 1000; % max_x value in x-axis in histogram to use data for bgnd peak fit.

disp(['xlim max for fit of bgnd peak = ' num2str(x_lim_for_fit)]) % print to command line.
fit_bgndPeak_x = binCentres2(binCentres2<x_lim_for_fit)'; % must be column vector for fitting.
fit_bgndPeak_y = freqCounts2(binCentres2<x_lim_for_fit)'; % must be column vector for fitting.
plot(fit_bgndPeak_x,fit_bgndPeak_y,'b-') % highlight in blue those values used later for Gaussian fit of bgnd peak.

% Fit part of the first peak to a Gaussian:
fun_Gauss = fittype('A*exp(-(x-x0)^2/(2*sigma^2))','independent','x'); % define Gaussian funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_A = 130; % set by eye from previous plot.
guess_sigma = 700; % set by eye from previous plot.
guess_x0 = 500; % set by eye from previous plot.
% Use coeffnames(fit_result) later to find out order of parameters.
options.StartPoint = [guess_A guess_sigma guess_x0]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [ ]; % Lower bounds for fit parameters. In order: A, sigma, x0.
options.Upper = [ ]; % Upper bounds for fit parameters. In order: A, sigma, x0.
[fit_result gof] = fit(fit_bgndPeak_x,fit_bgndPeak_y,fun_Gauss,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
A_fit = fit_param_values(1); 
sigma_fit = fit_param_values(2);
x0_fit = fit_param_values(3);
rsq_fit = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_A = errorSTDEV(1);
stDev_sigma = errorSTDEV(2);
stDev_x0 = errorSTDEV(3);

disp('Fit of background peak (on the right only up to its FWHM) to a Gaussian:')
disp(['A = ' num2str(A_fit) ' +- ' num2str(stDev_A)])
disp(['x0 = ' num2str(x0_fit) ' +- ' num2str(stDev_x0)])
disp(['sigma = ' num2str(sigma_fit) ' +- ' num2str(stDev_sigma)])
disp(['rsq_fit = ' num2str(rsq_fit)])

% Plot fit:
plot(fit_result,'k')
xlabel('Intensity')
ylabel('frequency')
hold off;


% Save also the fitted Gaussian to the first peak to the excel file:
x_fromFit = binCentres2';
y_fromFit = A_fit*exp(-(x_fromFit-x0_fit).^2./(2*sigma_fit^2));

xlswrite('mCherry_invitroResults',[x_fromFit,y_fromFit],'GaussFitFirstPeakInHistogram');

% Final result for in vitro for mCherry October 2012 is 500+-600.

% ---------------------
% % Mean value for the distribution data:
% mean_signal2 = sum(binCentres2.*freqCounts2)/sum(freqCounts2); % only positive y-values, all x-values in histogram.
% % Standard deviation value for the distribution data:
% stdev2 = sqrt(sum(((binCentres2-mean_signal2).^2).*freqCounts2)/sum(freqCounts2)); 
% 
% disp(' ') % empty line
% disp(['mean = ' num2str(mean_signal2)])
% disp(['std = ' num2str(stdev2)])

% ----------------------------------
% ----------------------------------
