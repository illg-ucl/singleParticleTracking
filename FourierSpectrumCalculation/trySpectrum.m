function trySpectrum()
%
% Created by Isabel Llorente-Garcia, Dic 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
%
% Function to try and find the best way to identify the spacing between
% periodic peaks in the histogram of intensity pair wise differences.
% Generate sum of periodic Gaussian peaks with noise and extract the
% spacing by calculating Fourier power spectrum.
%
%

%% Generate synthetic data for histogram of intensity pair-wise differences:
% with periodic peaks modelled as Gaussians.

x = 0:0.1:15; 

% Noise amplitude and vector:
noise_ampli = 20; % 0.5
noise_vector = noise_ampli*rand(1,length(x));
% Generate sum of random noise plus 5 gaussians:
% Vector of Gaussian centres:
period = 3;
Gauss_centres = [1 1+period 1+2*period 1+3*period 1+4*period];
% Vector of Gaussian widths:
Gauss_widths = 0.5*ones(1,5);
% Vector of Gaussian amplitudes:
Gauss_amplis = [12 27 18 0 100]; % [1 0.8 0.5 0.6 0.3]
% Gaussian functions, vectors:
gauss1 = Gauss_amplis(1)*gaussmf(x,[Gauss_widths(1) Gauss_centres(1)]);
gauss2 = Gauss_amplis(2)*gaussmf(x,[Gauss_widths(2) Gauss_centres(2)]);
gauss3 = Gauss_amplis(3)*gaussmf(x,[Gauss_widths(3) Gauss_centres(3)]);
gauss4 = Gauss_amplis(4)*gaussmf(x,[Gauss_widths(4) Gauss_centres(4)]);
gauss5 = Gauss_amplis(5)*gaussmf(x,[Gauss_widths(5) Gauss_centres(5)]);
% Add everything up:
y = noise_vector + gauss1 + gauss2 + gauss3 + gauss4 + gauss5;

% Create figure. 'position' vector is [left, bottom, width, height]. 
figure('color','white','units','inches','position',[0.2 1.5 18 10]); 
subplot(2,2,1);

plot(x,y)
title('Histogram of pair-wise intensity differences')
xlabel('Intensity pair-wise differences');
ylabel('probability');


%% Calculate spectrum (Fourier transform) of pair-wise difference distribution (of the histogram):

L = length(y); % length of probability y vector.
nFFT = 2^nextpow2(L); % number of points wanted in the Fourier Transform (it will either pad with zeros or truncate).
% Next power of 2 from length of vector (make the Fourier transform have a number of points which is a power of 2).
spectrum = fft(y,nFFT)/L; % fast fourier transform, spectrum of pair-wise difference probability distribution.

bin_separation = (x(end)-x(1))/(length(x)-1); % here the bin centres would be the x vector.
Fs = 1/bin_separation; % equivalent to sampling frequency, i.e. inverse of bin separation.
freq_axis = Fs/2*linspace(0,1,nFFT/2+1); % single-sided frequency axis.
% y = linspace(a,b,n) generates a row vector y of n points linearly spaced
% between and including a and b.
Istep_axis = 1./freq_axis; % axis with the intensity steps, i.e. the periodicity in terms of intensity jumps (period = 1/freq).

% Plot single-sided power spectrum:
subplot(2,2,2);
% plot the power spectrum in terms of the intensity_step:
power_spectrum = (2*abs(spectrum(1:nFFT/2+1))).^2; % single-sided power spectrum.
plot(Istep_axis,power_spectrum,'b-','LineWidth',1.5); 
xlim([0 15]);

title('Power Spectrum of pair-wise intensity difference histogram')
xlabel('Intensity period step')
ylabel('|Fourier Transform|^2')


%% Find peaks in the power spectrum of positive steps, to get the intensity step size:

[peak_ampli,peak_location] = findpeaks(power_spectrum); % see Matlab function 'findpeaks'.
Istep_peaks = Istep_axis(peak_location); % position of found peaks on the intensity_step axis.

% Put all peaks in a matrix, sort their amplitudes and show only the
% hightest 8:
peaks_to_print = [Istep_peaks' power_spectrum(peak_location)'];
sorted_peaks_to_print = sortrows(peaks_to_print,-2); % sort matrix in decreasing(-) order of 3rd column (power in spectrum).
selected_peaks_to_print = sorted_peaks_to_print(1:8,:); % highest 8 peaks.

% Mark found peaks on previous graph:
for i=1:length(selected_peaks_to_print)
    label_to_print = strcat('\leftarrow',num2str(selected_peaks_to_print(i,1),5));
    text(selected_peaks_to_print(i,1),selected_peaks_to_print(i,2),label_to_print,'FontSize',10);   
end


%% Second method: AUTOCORRELATION: 
% Try autocorrelation function to find periodicity:
% Autocorrelation as mathematical tool for finding repeating patterns,
% eg, periodic signal buried in noise, finding fundamental freq.
[autocor,x_lags] = xcorr(y,'unbiased'); % autocorrelation with 'unbiased' normalisation option (highlights peaks).
N = length(y);
% length of result autocor is twice length of y minus 1, i.e., 2*N-1.
subplot(2,2,3)
% calibrate autocorrelation axis:
x_axis_cor = x_lags*bin_separation;
% select only one half (other half is symmetrical) of the autocorrelation:
x_axis_autocor = x_axis_cor(N:2*N-1);
y_axis_autocor = autocor(N:2*N-1);
plot(x_axis_autocor,y_axis_autocor);
title('Autocorrelation method to find periodicity:');
xlabel('displacement (units of Intensity difference)')
ylabel('Autocorrelation function')

%% Find peaks in autocorrelation function:
min_dist_between_peaks = round((period-1)*length(x_axis_autocor)/max(x_axis_autocor));
[peak_ampli,peak_pos] = findpeaks(y_axis_autocor,'minpeakdistance',min_dist_between_peaks); % see Matlab function 'findpeaks'.
x_axis_peaks = x_axis_autocor(peak_pos); % position of found peaks on the intensity_step axis.

% Mark found peaks on previous graph:
for i=1:length(peak_pos)
    label_to_print = strcat('\leftarrow',num2str(x_axis_peaks(i),5));
    text(x_axis_peaks(i),peak_ampli(i),label_to_print,'FontSize',10);   
end


%% Calculate power spectrum (Fourier transform) of the autocorrelation:

% Use only one side of the autocorrelation:
% L2 = length(y_axis_autocor); % length of y vector.
% nFFT2 = 2^nextpow2(L2); % number of points wanted in the Fourier Transform (it will either pad with zeros or truncate).
% % Next power of 2 from length of vector (make the Fourier transform have a number of points which is a power of 2).
% spectrum2 = fft(y_axis_autocor,nFFT2)/L2; % fast fourier transform, spectrum of pair-wise difference probability distribution.
% 
% bin_separation2 = (x_axis_autocor(end)-x_axis_autocor(1))/(length(x_axis_autocor)-1); % here the bin centres would be the x vector.
% Fs2 = 1/bin_separation2; % equivalent to sampling frequency, i.e. inverse of bin separation.
% freq_axis2 = Fs2/2*linspace(0,1,nFFT2/2+1); % single-sided frequency axis.
% % y = linspace(a,b,n) generates a row vector y of n points linearly spaced
% % between and including a and b.
% Istep_axis2 = 1./freq_axis2; % axis with the intensity steps, i.e. the periodicity in terms of intensity jumps (period = 1/freq).
% 
% % Plot single-sided power spectrum:
% subplot(2,2,4);
% % plot the power spectrum in terms of the intensity_step:
% power_spectrum2 = (2*abs(spectrum2(1:nFFT2/2+1))).^2; % single-sided power spectrum.
% plot(Istep_axis2,power_spectrum2,'b-','LineWidth',1.5); 
% xlim([0 15]);
% 
% title('Power Spectrum of autocorrelation')
% xlabel('Intensity period step')
% ylabel('|Fourier Transform|^2')

% Use both sides of the autocorrelation:
L2 = length(autocor); % length of y vector.
nFFT2 = 2^nextpow2(L2); % number of points wanted in the Fourier Transform (it will either pad with zeros or truncate).
% Next power of 2 from length of vector (make the Fourier transform have a number of points which is a power of 2).
spectrum2 = fft(autocor,nFFT2)/L2; % fast fourier transform, spectrum of pair-wise difference probability distribution.

bin_separation2 = (x_axis_cor(end)-x_axis_cor(1))/(length(x_axis_cor)-1); % here the bin centres would be the x vector.
Fs2 = 1/bin_separation2; % equivalent to sampling frequency, i.e. inverse of bin separation.
freq_axis2 = Fs2/2*linspace(0,1,nFFT2/2+1); % single-sided frequency axis.
% y = linspace(a,b,n) generates a row vector y of n points linearly spaced
% between and including a and b.
Istep_axis2 = 1./freq_axis2; % axis with the intensity steps, i.e. the periodicity in terms of intensity jumps (period = 1/freq).

% Plot single-sided power spectrum:
subplot(2,2,4);
% plot the power spectrum in terms of the intensity_step:
power_spectrum2 = (2*abs(spectrum2(1:nFFT2/2+1))).^2; % single-sided power spectrum.
plot(Istep_axis2,power_spectrum2,'b-','LineWidth',1.5); 
xlim([0 15]);

title('Power Spectrum of autocorrelation')
xlabel('Intensity period step')
ylabel('|Fourier Transform|^2')


%% Find peaks in the power spectrum of the autocorrelation function:

[peak_ampli2,peak_location2] = findpeaks(power_spectrum2); % see Matlab function 'findpeaks'.
Istep_peaks2 = Istep_axis2(peak_location2); % position of found peaks on the intensity_step axis.

% Put all peaks in a matrix, sort their amplitudes and show only the
% hightest 8:
peaks_to_print2 = [Istep_peaks2' power_spectrum2(peak_location2)'];
sorted_peaks_to_print2 = sortrows(peaks_to_print2,-2); % sort matrix in decreasing(-) order of 2ND column (power in spectrum).
selected_peaks_to_print2 = sorted_peaks_to_print2(1:8,:); % highest 8 peaks.

% Mark found peaks on previous graph:
for i=1:length(selected_peaks_to_print2)
    label_to_print2 = strcat('\leftarrow',num2str(selected_peaks_to_print2(i,1),5));
    text(selected_peaks_to_print2(i,1),selected_peaks_to_print2(i,2),label_to_print2,'FontSize',10);   
end
