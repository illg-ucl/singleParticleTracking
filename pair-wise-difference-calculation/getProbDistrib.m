function [binCentres freqCounts] = getProbDistrib(column_vector,points_per_bin,display_figures)
%
% Created by Isabel Llorente Garcia. March 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Calculation of probability distribution (or histogram) from an input
% vector and optional display of their histograms.
%
% Inputs: 
% column_vector: needs to be a COLUMN vector. Contains values we want to
% bin and get the probability of, the vector we want to make a histogram of. 
% (typically contains all pair-wise differences).
% points_per_bin: related to the bin size for the histograms/or
% probability distribution. 
% (Default is 20 when there is a large unsubtracted background and around 2 when the signal has
% is background-subtracted.)
% display_figures: 1 if you want to show plots, 0 if not.
%
% Outputs:
% binCentres: row vector containing the resulting bin centres of the full
% histogram of all values in the column vector.
% freqCounts: row vector containing the occurrences (or frequency) for each
% value in the column_vector.
% The first plot shows the full histogram.
% The second plot shows the highest end of the full histogram.
% The third plot shows the lower end of the full histogram.


% Calculate and show pairwise difference probability distribution (plot histograms):
nbins = round(length(column_vector)/points_per_bin); % chosen number of bins which seems to work fine in most cases. Cannot be too small!
[freqCounts,binCentres] = hist(column_vector,nbins); % histogram. 

% Display figures only if input display_figures is 1:
if display_figures == 1 
    
    figure('position',[1200 100 500 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
    subplot(3,1,1) % full histogram.
    bar(binCentres,freqCounts,'r'); % plot a bar graph of the full histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    ylim([0 1.1*max(freqCounts)]); % re-scale vertical axis.
    
    % Separate components in previous histogram by changing range to plot.
    % If the signal is not background-subtracted: the lower end of the large 
    % positive values of pair-wise differences in the histogram
    % corresponds to the background level (if no background subtraction has
    % been applied).
    % The smaller positive and negative values of pair-wise differences in the
    % histogram correspond to either photobleaching jumps, blinking or noise.
    subplot(3,1,2) % Background component:
    bar(binCentres,freqCounts,'b'); % plot same full histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title('background component');
    xlim([0.9*max(binCentres) 1.1*max(binCentres)]); % re-scale horizontal axis.
    
    subplot(3,1,3) % Photobleaching/blinking/noise component:
    bar(binCentres,freqCounts,'g'); % plot same full histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title('photobleaching/blinking/noise component');
    xlim([1.1*min(binCentres) -1.1*min(binCentres)]); % re-scale horizontal axis.
    
end