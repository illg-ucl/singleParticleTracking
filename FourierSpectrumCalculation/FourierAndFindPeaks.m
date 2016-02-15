function [power_spectrum_x power_spectrum_y spectrum_peaks_x spectrum_peaks_y] = FourierAndFindPeaks(binCentres,freqCounts,display_figure,Xmax_inplot)
% 
% Created by Isabel Llorente Garcia. March 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Calculate spectrum (Fourier transform) of histogram data (eg, of pair-wise difference
% distribution) and find peaks in the power spectrum 
% (to get the intensity step size).
% 
% Inputs: binCentres and freqCounts are usually the outputs of function
% fullPwD.m, they are the x and y data for the histogram of full pair-wise
% differences of intensity.
% display_figure: 1 for yes, 0 for no.
% Xmax_inplot: maximum value of x-axis (Istep value) to plot.
%
% Outputs: all column vectors:
% power_spectrum_x: x-axis of power spectrum (Intensity step values).
% power_spectrum_y: y-axis of power spectrum (spectral power).
% spectrum_peaks_x: x-axis (Istep value) position of peaks found in spectrum,
% sorted in terms of increasing spectral power
% spectrum_peaks_y: y-axis (spectral power) position of peaks found in spectrum,
% sorted in terms of increasing spectral power.


%% Calculate spectrum (Fourier transform) of pair-wise difference distribution:

L = length(freqCounts); % length of probability (freq counts) for pwd.
nFFT = 2^nextpow2(L); % number of points wanted in the Fourier Transform (it will either pad with zeros or truncate).
% Next power of 2 from length of vector (make the Fourier transform have a number of points which is a power of 2).
spectrum_pwd_distrib = fft(freqCounts,nFFT)/L; % fast fourier transform, spectrum of pair-wise difference probability distribution.

bin_separation = (binCentres(end)-binCentres(1))/(length(binCentres)-1);
Fs = 1/bin_separation; % equivalent to sampling frequency, i.e. inverse of bin separation.
freq_axis = Fs/2*linspace(0,1,nFFT/2+1); % single-sided frequency axis in Hz.
% y = linspace(a,b,n) generates a row vector y of n points linearly spaced
% between and including a and b.

warning('off','MATLAB:divideByZero'); % this turns off warning if Divide By Zero for the following. 
% Use warning('query','last') on command window to find warning identifier.
Istep_axis = 1./freq_axis; % axis with the intensity steps, i.e. the periodicity in terms of intensity jumps (period = 1/freq).
% If first point in freq_axis is zero, first point in Istep_axis is Inf,
% but it's ok.

% Calculate the power spectrum:
power_spectrum = (2*abs(spectrum_pwd_distrib(1:nFFT/2+1))).^2; % single-sided power spectrum.

% % Error control needed for Matlab version R2011b (otherwise, findpeaks(0) gives an error):
% if power_spectrum == 0
%     power_spectrum = [0 0 0 0];
% end


%% Find peaks in the power spectrum of positive steps, to get the intensity step size:

[peak_ampli,peak_location] = findpeaks(power_spectrum); % see Matlab function 'findpeaks'.
warning('off','signal:findpeaks:noPeaks'); % this turns off warning if no peaks found. 
% Use warning('query','last') on command window to find warning identifier.

if ~isempty(peak_location)  % if peaks are found (at least one)
    
    Istep_peaks = Istep_axis(peak_location); % position of found peaks on the intensity_step axis.
    
    % Put all peaks in a matrix, sort their amplitudes and show only the
    % hightest 8:
    peaks_to_print = [Istep_peaks' power_spectrum(peak_location)'];
    sorted_peaks_to_print = sortrows(peaks_to_print,-2); % sort matrix in decreasing(-) order of 2nd column (power in spectrum).
    numOfPeaksToShow = min(size(sorted_peaks_to_print,1),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    selected_peaks_to_print = sorted_peaks_to_print(1:numOfPeaksToShow,:); % select highest peaks.
    
    % Show results graphically:
    if display_figure==1
        % Plot single-sided power spectrum versus intensity_step values:
        figure
        % plot(freq_axis,power_spectrum,'b-','LineWidth',2.5);
        plot(Istep_axis,power_spectrum,'b-','LineWidth',1.5);
        xlim([0 Xmax_inplot]);
        % ylim([0 1]);
        % change horizontal label tick marks:
        % freq_ticks = get(gca,'XTick'); % regularly spaced ticks which would appear for the previous frequency axis.
        % Istep_ticks = 1./freq_ticks; % intensity_step period axis ticks. Now need to transform this into a cell array of strings:
        % Istep_tick_labels = cell(1,length(Istep_ticks)); % initialise empty cell of strings for I_step_period tick labels.
        % for j =1:length(Istep_ticks)
        %     Istep_tick_labels{j} = num2str(Istep_ticks(j),5); % the ,5 gives two digits of precission in the num2str conversion.
        % end
        % set(gca,'XTickLabel',Istep_tick_labels)
        title('Power Spectrum of pair-wise intensity difference probability distribution')
        xlabel('Intensity period step')
        ylabel('|Fourier Transform|^2')
        % Mark found peaks on previous graph:
        for i=1:size(selected_peaks_to_print,1)
            label_to_print = strcat('\leftarrow',num2str(selected_peaks_to_print(i,1),3)); % print 3 significant figures.
            text(selected_peaks_to_print(i,1),selected_peaks_to_print(i,2),label_to_print,'FontSize',10);
        end
    end
    
    % Outputs:
    power_spectrum_x = Istep_axis'; % column vector, x-axis of power spectrum (Intensity step values).
    power_spectrum_y = power_spectrum'; % column vector, y-axis of power spectrum (spectral power values).
    spectrum_peaks_x = sorted_peaks_to_print(:,1); % column vector with intensity steps, i.e., found peaks in spectrum.
    spectrum_peaks_y = sorted_peaks_to_print(:,2); % power (y-axis) of found peaks in spectrum.
else
    % Outputs:
    power_spectrum_x = Istep_axis'; % column vector, x-axis of power spectrum (Intensity step values).
    power_spectrum_y = power_spectrum'; % column vector, y-axis of power spectrum (spectral power values).
    spectrum_peaks_x = [];
    spectrum_peaks_y = [];
end

