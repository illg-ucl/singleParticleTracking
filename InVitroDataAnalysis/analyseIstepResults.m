function analyseIstepResults(Isteps_vector,nbins1,limit1,nbins2,limit2,method)
% 
% Created by Isabel Llorente-Garcia. April 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Function to analyse (in two ways) the obtained results for values of
% intensity steps (from in vitro data of fluorophores, for instance).
%
% Inputs: 
% "Isteps" is a column vector containing all the intensity-step values.
% "nbins1": number of bins in first histogram of Isteps. (100).
% "limit1": maximum value (intensity units) in x-axis of spectrum of histogram of Isteps.
% "nbins2": number of bins in second histogram of intensity pair-wise differences. (100).
% "limit2": maximum value (intensity units) in x-axis of spectrum of histogram of pair-wise differences of Isteps.
%

% disp(['Number of spots analysed: ',num2str(length(Isteps_vector))]);
disp(['Mean Istep: ',num2str(mean(Isteps_vector))]);
disp(['Standard deviation Istep: ',num2str(std(Isteps_vector))]);
disp(['Median Istep: ',num2str(median(Isteps_vector))]);
% disp(['Mode Istep: ',num2str(mode(Isteps_vector))]);

%% METHOD 1:

if method == 1
    % Calculate and plot histogram of results:
    [freqCounts_0,binCentres_0] = hist(Isteps_vector,nbins1); % good number of bins is 100.
    figure; bar(binCentres_0,freqCounts_0,'g'); % plot a bar graph of the full histogram.
    % xlabel('Intensity-steps from spots');
    xlabel('Intensities from all spots');
    ylabel('frequency');
    title('Histogram of Intensities');

    
%     % Find the parameters of the gamma distribution which best fits the data:
%     [params_gamma, pci] = gamfit(Isteps_vector);
%     % For a gamma distribution with parameters a and b, gampdf(x,a,b), params_gamma(1)
%     % gives "a" and params_gamma(2) gives "b". "pci" gives the 95% confidence intervals
%     % for "a" and "b" in each column, respectively. The 95% confidence intervals for "a" is the first column, pci(1,1) to
%     % pci(2,1), and for "b" it is the second column: pci(1,2) to pci(2,2).
%     [mean_gamma,variance_gamma] = gamstat(params_gamma(1),params_gamma(2)); % returns mean of and variance for the gamma distribution with parameters specified by params_gamma(1) and params_gamma(2).
%     disp(['params gamma distrib:  ',num2str(params_gamma)]);
%     disp(['mean_gamma: ',num2str(mean_gamma)]);
%     disp(['std_gamma: ',num2str(sqrt(variance_gamma))]);
%     disp(['mode gamma: ',num2str((params_gamma(1)-1)*params_gamma(2))]);
% 
%     % Gamma probability density distribution function:
%     gamma_distrib_pdf = gampdf(binCentres_0,params_gamma(1),params_gamma(2));
%     
%     bin_size = binCentres_0(2)-binCentres_0(1);
%     normalisation_const = sum(bin_size*freqCounts_0); % normalisation constant, total area under prob density distrib curve.
%     freqCounts_normalised = freqCounts_0/normalisation_const;
%     figure; 
%     bar(binCentres_0,freqCounts_normalised,'g'); % plot a bar graph of the full histogram.
%     hold on;
%     plot(binCentres_0,gamma_distrib_pdf,'-k');
%     title('normalised prob distrib and fit to gamma distrib');
%     hold off;
%     
%     % Fit histogram (not normalised) of intensity levels to a gamma distribution:
%     % (because if we cut the distribution by having a threshold>0, the previous
%     % command gamfit does not work well to try and fit the distribution,
%     % due to wrong normalisation.)
%     fun_to_fit = fittype('gampdf(x,a,b)*c','independent','x'); % define exponential funtion to fit to, with 'x' as independent variable;
%     options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
%     % Guesses for fit parameters:
%     options.StartPoint = [params_gamma(1) params_gamma(2) normalisation_const]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
%     [fit_result gof] = fit(binCentres_0',freqCounts_0',fun_to_fit,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
%     % fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
%     fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
%     a_fit = fit_param_values(1);
%     b_fit = fit_param_values(2);
%     c_fit = fit_param_values(3);
%     rsq_fit = gof.rsquare; % rsquare coefficient of fit.
%     errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
%     errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
%     stDev_a = errorSTDEV(1);
%     stDev_b = errorSTDEV(2);
%     stDev_c = errorSTDEV(3);
%     
%     disp(' ') % empty line
%     disp('Fit of intensity prob density to gamma distrib: ')
%     disp([' a = ',num2str(a_fit),' +- ',num2str(stDev_a),';   b = ',num2str(b_fit),' +- ',num2str(stDev_b)])
%     disp([' c = ',num2str(c_fit),' +- ',num2str(stDev_c),';    rsqr = ',num2str(rsq_fit)])
%     
%     figure;
%     bar(binCentres_0,freqCounts_0,'r'); % plot a bar graph of the full histogram.
%     hold on;
%     plot(fit_result,'-k');
%     xlabel('Intensities from all spots');
%     ylabel('frequency');
%     title('Prob density and fit to gamma distrib');
%     hold off;
%     
%     disp('Gamma stats from new fit: ')
%     [mean_gamma_new,variance_gamma_new] = gamstat(a_fit,b_fit); % returns mean of and variance for the gamma distribution with parameters specified by params_gamma(1) and params_gamma(2).
%     disp(['new params gamma distrib:  ',num2str([a_fit b_fit])]);
%     disp(['mean_gamma_new: ',num2str(mean_gamma_new)]);
%     disp(['std_gamma_new: ',num2str(sqrt(variance_gamma_new))]);
%     disp(['mode gamma_new: ',num2str((a_fit-1)*b_fit)]);
    
    
    % Fourier transform, spectrum of previous histogram:
    [ps_x_0 ps_y_0 ps_peaks_x_0 ps_peaks_y_0] = FourierAndFindPeaks(binCentres_0,freqCounts_0,1,limit1);
    title('Spectrum of histogram of Intensity steps');
    
end


%% METHOD 2:

if method == 2
    % Better method (gives clearer peak), use this method:
    % calculate pairwise differences first, then histogram of them:
    [binCentres freqCounts] = fullPwD(Isteps_vector,nbins2,0);
    % one but last param is number of bins (100).
    figure; bar(binCentres,freqCounts,'r');
    xlabel('Pair-wise differences of intensity-steps');
    ylabel('frequency');
    title('Histogram of pair-wise differences of Intensity steps');
    % Fourier transform, spectrum of previous histogram:
    [ps_x ps_y ps_peaks_x ps_peaks_y] = FourierAndFindPeaks(binCentres,freqCounts,1,limit2);
    title('Spectrum of histogram of pair-wise differences of Intensity steps');
    % last two params are display_figure and xmax to plot spectrum.
end