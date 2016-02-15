
% runToAnalyse.m
% In vitro data analysis. GFP
% --------------------------------
%
% limitsROI = [250 500 250 500];
% limitsYzoom = [-2000 6000];

% params.SNR_min = 2; 
% params.rsq_min = 0.2;

% r554 = analyseFixedBrightSpots('554',15,100,'');
% % for params.x_limit_spectr = 2000;
% 
% r556 = analyseFixedBrightSpots('556',18,100,'');
% % for params.x_limit_spectr = 2000;
% 
% r557 = analyseFixedBrightSpots('557',22,100,'');
% % for params.x_limit_spectr = 2000;
% 
% r558 = analyseFixedBrightSpots('558',11,100,'');
% % for params.x_limit_spectr = 2000;
% 
% 
% % % r554bis = analyseFixedBrightSpots('554',15,100,'bis');
% % % % for params.x_limit_spectr = 3000;
% % % 
% % % r556bis = analyseFixedBrightSpots('556',18,100,'bis');
% % % % for params.x_limit_spectr = 3000;
% % % 
% % % r557bis = analyseFixedBrightSpots('557',22,100,'bis');
% % % % for params.x_limit_spectr = 3000;
% % % 
% % % r558bis = analyseFixedBrightSpots('558',11,100,'bis');
% % % % for params.x_limit_spectr = 3000;
% 
% % Save results:
% save 'results' 'r*' % save all results in a .mat file in current folder.



% ------------------------
%% Put together intensities of all spots (or selected spots only) and make histogram of them, with threshold:

load results.mat;

disp(' '); % empty line
disp('For ALL spots:');
allImages = {r554,r556,r557,r558}; % results from all images.
intensities_vector_allImages = []; % initialise vector to accummulate intensity values from all images;
threshold = -Inf; % only accummulate intensity values above this threshold (use -Inf to use all values).   
for j =1:length(allImages)
    
    r_toAnalyse = allImages{j};
    
    r_Ivectors = [r_toAnalyse{1}.integrated_Ispot_bgndSubtracted]; % matrix in which each column is the intensity of a given spot, each row is frame no.
    intensities_vector = []; % initialise vector to accummulate intensity values;
    for i = 1:size(r_Ivectors,2)
        % go through selected spots in these image:
        intensity_values = r_Ivectors(:,i); % column vector.
        intensities_vector = [intensities_vector; intensity_values(intensity_values > threshold)]; % append values larger than threshold;
    end
    
    intensities_vector_allImages = [intensities_vector_allImages; intensities_vector];
    
end
% % Check that results are robust to change in the number of histogram bins, as they should:
analyseIstepResults(intensities_vector_allImages,250,1400,250,1400,1)
analyseIstepResults(intensities_vector_allImages,500,1400,500,1400,1)
analyseIstepResults(intensities_vector_allImages,1000,1400,1000,1400,1)

% % Check also all the spots corresponding to a single image sequence (vector intensities_vector) (put a breakpoint inside loop):
% % Calculate and plot histogram of all intensity values and also of their pair-wise differences:
% analyseIstepResults(intensities_vector,100,4000,100,4000,1) % make the x-axis max limit of the Fourier spectra equal to 4000.
% analyseIstepResults(intensities_vector,100,2000,100,2000,1) % make the x-axis max limit of the Fourier spectra equal to 2000 (zoom in).
% analyseIstepResults(intensities_vector,100,1000,100,1000,1) % make the x-axis max limit of the Fourier spectra equal to 1000 (zoom further in).


% Save to existing output excel file: 
% CAREFUL!!! It will overwrite existing file!
xlswrite('invitroResults',intensities_vector_allImages,'all_I_values'); % save to sheet 'all_I_values_noThreshold'.
% Calculate and plot histogram of all values:
% x and y values of the histogram plot:
[freqCounts,binCentres] = hist(intensities_vector_allImages,250);
[freqCounts1,binCentres1] = hist(intensities_vector_allImages,500);
[freqCounts2,binCentres2] = hist(intensities_vector_allImages,1000);
% Save to existing output excel file: 
% CAREFUL!!! It will overwrite existing file!
xlswrite('invitroResults',[binCentres',freqCounts'],'histogram x-y all_250bars');
xlswrite('invitroResults',[binCentres1',freqCounts1'],'histogram x-y all_500bars');
xlswrite('invitroResults',[binCentres2',freqCounts2'],'histogram x-y all_1000bars');

% Plot histogram in green:
bar(binCentres2,freqCounts2,'g');
hold on;

% Take only the left hand side of main peak (highlighted in blue on plot) and only up to approx. the FWHM to its right hand side and fit it to a gaussian.
% We assume the distribution is Gaussian. Fit shown as a black curve.
x_lim_for_fit = 0; % max_x value in x-axis in histogram to use data for bgnd peak fit.

disp(['xlim max for fit of bgnd peak = ' num2str(x_lim_for_fit)]) % print to command line.
fit_bgndPeak_x = binCentres2(binCentres2<x_lim_for_fit)'; % must be column vector for fitting.
fit_bgndPeak_y = freqCounts2(binCentres2<x_lim_for_fit)'; % must be column vector for fitting.
plot(fit_bgndPeak_x,fit_bgndPeak_y,'b-') % highlight in blue those values used later for Gaussian fit of bgnd peak.

% Fit part of the first peak to a Gaussian:
fun_Gauss_0 = fittype('A*exp(-(x-0)^2/(2*sigma^2))','independent','x'); % define Gaussian funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_A = 30; % set by eye from previous plot.
guess_sigma = 200; % set by eye from previous plot.
% guess_x0 = 0; % set by eye from previous plot.
% Use coeffnames(fit_result) later to find out order of parameters.
options.StartPoint = [guess_A guess_sigma]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [ ]; % Lower bounds for fit parameters. In order: A, sigma, x0.
options.Upper = [ ]; % Upper bounds for fit parameters. In order: A, sigma, x0.
[fit_result gof] = fit(fit_bgndPeak_x,fit_bgndPeak_y,fun_Gauss_0,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
A_fit = fit_param_values(1); 
sigma_fit = fit_param_values(2);
% x0_fit = fit_param_values(3);
rsq_fit = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_A = errorSTDEV(1);
stDev_sigma = errorSTDEV(2);
% stDev_x0 = errorSTDEV(3);

disp('Fit of background peak (on the right only up to its FWHM) to a Gaussian:')
disp(['A = ' num2str(A_fit) ' +- ' num2str(stDev_A)])
% disp(['x0 = ' num2str(x0_fit) ' +- ' num2str(stDev_x0)])
disp(['sigma = ' num2str(sigma_fit) ' +- ' num2str(stDev_sigma)])
disp(['rsq_fit = ' num2str(rsq_fit)])

% Plot fit:
plot(fit_result,'k')
xlabel('Intensity')
ylabel('frequency')
hold off;


% Save also the fitted Gaussian to the first peak to the excel file:
x_fromFit = binCentres2';
y_fromFit = A_fit*exp(-(x_fromFit-0).^2./(2*sigma_fit^2));

xlswrite('invitroResults',[x_fromFit,y_fromFit],'GaussFitFirstPeakInHistogram');





% Subtract the fitted peak around zero and plot what remains in red:
figure;
y_afterRemoveZeroPeak = freqCounts2-y_fromFit';
plot(binCentres2,y_afterRemoveZeroPeak,'r');
xlabel('Intensity')
ylabel('frequency')
title('After subtracting peak centred on zero')
hold on;

% Save to excel file:
xlswrite('invitroResults',[x_fromFit,y_afterRemoveZeroPeak'],'AfterSubtractFirstPeak');


% Fit the remaining histogram after subtracting peak around zero to a Gaussian:
% Fit only the positive y-values:
fit_Peak_x = binCentres2(y_afterRemoveZeroPeak>0)'; % must be column vector for fitting.
fit_Peak_y = y_afterRemoveZeroPeak(y_afterRemoveZeroPeak>0)'; % must be column vector for fitting.

fun_Gauss = fittype('A*exp(-(x-x0)^2/(2*sigma^2))','independent','x'); % define Gaussian funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_A = 30; % set by eye from previous plot.
guess_sigma = 300; % set by eye from previous plot.
guess_x0 = 400; % set by eye from previous plot.
% Use coeffnames(fit_result) later to find out order of parameters.
options.StartPoint = [guess_A guess_sigma guess_x0]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [ ]; % Lower bounds for fit parameters. In order: A, sigma, x0.
options.Upper = [ ]; % Upper bounds for fit parameters. In order: A, sigma, x0.
[fit_result gof] = fit(fit_Peak_x,fit_Peak_y,fun_Gauss,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
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

disp('Fit remaining histogram (after subtracting peak centred on zero) to a Gaussian:')
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
y_fromFit_afterSubtract = A_fit*exp(-(x_fromFit-x0_fit).^2./(2*sigma_fit^2));

xlswrite('invitroResults',[x_fromFit,y_fromFit_afterSubtract],'GaussFitAfterSubtractFirstPeak');


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
