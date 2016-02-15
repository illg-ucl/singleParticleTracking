function rescaleI0probDensityDistrib()
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Compare and re-scale probability density distributions if initial intensities (I0, prop to stoichiometry) for different data
% sets.
% Import the results of all values of I0 from the analysis of many tracks:

% ------------------------ 
% Import data:
% For the single label
path1 = 'C:\Isabel\ExperimData\HeikoData\ATPase-GFP\ATPase-GFP TIRF\good\ATPase-GFP_InitialI_bottom2p8\results_allTracks\results_ATPase-GFP_bottom.mat';
results1 = load(path1); % This is a structure with fields:
%               results: [1x... struct]
%     allTrack_results1: [1x1 struct]
%     allTrack_results2: [1x1 struct]
%       allTrack_inputs: [1x1 struct]
x1 = results1.allTrack_results2.I0_axis;
y1 = results1.allTrack_results2.I0_prob_density;

% For the dual label:
path2 = 'C:\Isabel\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\cybD-mCherry-ATPase-GFp_InitialI_bottom4p4\results_allTracks\results_cybD-mCherry-ATPase-GFp_bottom.mat';
results2 = load(path2); % This is same type of structure as before.
x2 = results2.allTrack_results2.I0_axis;
y2 = results2.allTrack_results2.I0_prob_density;

% ------------------------ 
% Plot both prob density distributions together.
plot(x1,y1,'c') % single label, CYAN.
hold on; 
plot(x2,y2,'b') % dual label, in BLUE.
hold off;
xlabel('Initial Intensity');
ylabel('Probability density');
legend('single label','dual label');

% Wait for users input (give time to look at plot):
input('Save figure, maybe? Press any key to continue...');

% ------------------------ 
% Re-scale one distribution (single label one) to fit the other one@

first_scaling_factor = 0.5;
last_scaling_factor = 1.2;
step_iteration = 0.01;
scaling_factors = (first_scaling_factor:step_iteration:last_scaling_factor);
n_iterations = length(scaling_factors);

sos = zeros(n_iterations,1); % initialise empty variable to store sum of squares for each iteration.

for i=1:n_iterations
    
    disp(['scaling_factor: ' num2str(scaling_factors(i))]);
    
    x1_rescaled = scaling_factors(i)*x1;
    y1_rescaled = y1/scaling_factors(i); % rescale y axis too to keep the whole area under the curve (whole prob) equal to 1.
    
    % Resample y1_rescaled so that its x-axis has the same values as x2:
    % Resampling does not affect the normalisation to a total probability of 1.
    % y1_rescaled_resampled = interp1q(x1_rescaled,y1_rescaled,x2);
    y1_rescaled_resampled = interp1(x1_rescaled,y1_rescaled,x2,'linear',0);  % use 0 for out of range values.
    
    % Now the new distribution is:
    x3 = x2;
    y3 = y1_rescaled_resampled;
    
%     % Graphic check:
%     figure; plot(x1_rescaled,y1_rescaled)
%     hold; plot(x3,y3,'m');
    
    % Compare the sum of squared differences for all points in distributions given by (x3,y3) and (x2,y2):
    sos(i) = sum((y2-y3).^2);
    
    % Plot the original distribution (x2,y2) and the re-scaled one (x3,y3):
    plot(x2,y2,'b') % dual label, BLUE.
    hold on;
    plot(x3,y3,'m') % single label, re-scaled, in MAGENTA.
    hold off;
    xlabel('Initial Intensity');
    ylabel('Probability density');
    pause(0.3); % pause for 1 second.
    
end

% Plot resulting sum of squares versus scaling factor:
figure;
plot(scaling_factors,sos,'-k'); 
xlabel('scaling factor');
ylabel('sum of squares');

% ---------------------------------
% Result:
% find scaling factor which minimises the sum of squares:
scaling_factor_result = scaling_factors(sos==min(sos)); 
% Print result to command line:
disp(' '); % empty line
disp(['scaling_factor_result: ' num2str(scaling_factor_result)]);
% Plot overlay of original distribution and final re-scaled one:
x1_rescaled = scaling_factor_result*x1;
y1_rescaled = y1/scaling_factor_result; 
figure;
plot(x2,y2,'b') % original dual label, BLUE.
hold on;
plot(x1_rescaled,y1_rescaled,'m') % re-scaled single label, in MAGENTA.
hold off;
xlabel('Initial Intensity');
ylabel('Probability density');
legend('dual label','re-scaled single label');

% Chi-square test on the match of the two distributions: ((observed minus expected)^2/expected)
y1_rescaled_resampled = interp1(x1_rescaled,y1_rescaled,x2,'linear',0); % re-sample first to have numbers at same x-values for both distribs.
chi_square_0 = sum(((y1_rescaled_resampled-y2).^2)./y2);
% Compare the chi_square_0 value obtained with that for a Chi-Square distrib
% of a number of degrees of freedom equal to the no. of points compared minus one (length(x2)-1), and at the prob_limit probability level:
prob_limit = 0.01; % 1% (can be 0.01 for a small prob acceptance).
chi_square_limit = chi2inv(1-prob_limit,length(x2)-1);
% The chi_square_0 value should be lower than the chi_square_limit value
% for the two distributions to agree to a probability level of prob_limit.
chi_square_50 = chi2inv(0.5,length(x2)-1); % 50 per cent prob limit
% Cummulative probability distrib function:
cumul_prob_0 = chi2cdf(chi_square_0,length(x2)-1);
% The value of cumul_prob_0 is the cumulative probability that a value in
% the corresponding Chi-square distribution is found below the obtained
% chi_square_0 value. This should be small enough compared to the chosen
% prob_limit value. Ideally, it should be close to 0.5.
% 1-cumul_prob_0 should be above the prob_limit.
disp(' ');
disp(['chi_square_0: ' num2str(chi_square_0)]);
disp(['chi_square for prob ' num2str(prob_limit) ': ' num2str(chi_square_limit)]); % Chi-square value for 50% cumulative probability.
disp(['chi_square for 50% prob ' num2str(chi_square_50)]); % Chi-square value for 50% cumulative probability.
disp(['cumul_prob_0: ' num2str(cumul_prob_0)]);
disp(['1-cumul_prob_0: ' num2str(1-cumul_prob_0)]);
disp(' ');
