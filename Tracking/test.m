n =1; % Initialise index for position in final result
for i = 1:size(Inormalised_all_ATPase_GFP,1) % loop through rows in matrix of intensity traces Inormalised_all:
    
    Inormalised_row0 = Inormalised_all_ATPase_GFP(i,:); % take that row.
    Inormalised_row = Inormalised_row0(Inormalised_row0 ~= 0); % exclude values equal to zero (empty data positions).
    if length(Inormalised_row)>4 % if there are at least 4 elements different from zero in that row:
        Inormalised_avg(n) = mean(Inormalised_row); % mean intensity value for that row, ie, for that time (relative time).
        Inormalised_error(n) = std(Inormalised_row)/sqrt(length(Inormalised_row)); % use standard deviation/sqrt(num points) as error of the mean.
        n = n+1;
    end
    
end

% [Inormalised_avg' Inormalised_error']   % display result
% % Results for later plot:
timeAxis = 0.04.*[0:1:length(Inormalised_avg)-1]; % row vector


%% Fit intensity data to an exponentially decaying function (with no offset):

% We fit the averaged and normalised integrated spot intensity (bgnd subtracted) to an exponential:
IforFit = Inormalised_avg'; % column vector
tforFit = timeAxis'; %  column vector,  time relative to start of track.
exp_no_offset = fittype('I0*exp(-t/tau)','independent','t'); % define exponential funtion to fit to, with 't' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = 1; % since I has been normalised.
guess_tau = 2; % in seconds;
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 0]; % Lower bounds for fit parameters. In order: I0, tau.
options.Upper = [2 30]; % Upper bounds for fit parameters. In order: I0, tau.
[fit_result_I gof] = fit(tforFit,IforFit,exp_no_offset,options); % fit_result_I contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% fit_param_names = coeffnames(fit_result_I); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result_I); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
I0_fit = fit_param_values(1); % I0 intensity value from fit.
tau_fit = fit_param_values(2); % tau from fit.
rsq_fit_I = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result_I,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_I0 = errorSTDEV(1);
stDev_tau = errorSTDEV(2);

disp(' ') % empty line
disp('Intensity vs time exponential fit (with no offset) result: ') 
disp([' I0 = ',num2str(I0_fit),' +- ',num2str(stDev_I0),';   tau = ',num2str(tau_fit),' +- ',num2str(stDev_tau),' s.',';   rsq = ',num2str(rsq_fit_I)]) 

% results from exponential fit with no offset of intensity vs time:
results_I_fits.I0_fit = I0_fit;
results_I_fits.stDev_I0 = stDev_I0;
results_I_fits.I0_fit_percentError = 100*stDev_I0/I0_fit;
results_I_fits.tau_fit = tau_fit;
results_I_fits.stDev_tau = stDev_tau;
results_I_fits.tau_fit_percentError = 100*stDev_tau/tau_fit;
results_I_fits.rsq_fit_I = rsq_fit_I;


%% Fit intensity data to an exponentially decaying function (with offset "_wo"):

exp_with_offset = fittype('I0*exp(-t/tau)+Ioffset','independent','t'); % define exponential funtion to fit to, with 't' as independent variable and 'tstart' as a fixed parameter (constant);
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = 1; 
guess_Ioffset = 0;
guess_tau = 2; 
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_Ioffset guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 -0.5 0]; % Lower bounds for fit parameters. In order: I0, Ioffset, tau.
options.Upper = [2 0.5 30]; % Upper bounds for fit parameters. In order: I0, Ioffset, tau.
try % error control in case fit fails
    [fit_result_I_wo gof] = fit(tforFit,IforFit,exp_with_offset,options); % fit_result_I_wo contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
    % fit_param_names = coeffnames(fit_result_I_wo); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
    fit_param_values = coeffvalues(fit_result_I_wo); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
    I0_fit_wo = fit_param_values(1); % I0 intensity value from fit.
    Ioffset_fit_wo = fit_param_values(2); % Ioffset from fit.
    tau_fit_wo = fit_param_values(3); % tau from fit.
    rsq_fit_I_wo = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_result_I_wo,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    stDev_I0_wo = errorSTDEV(1);
    stDev_Ioffset_wo = errorSTDEV(2);
    stDev_tau_wo = errorSTDEV(3);
catch ME1
    fit_result_I_wo = [0];
    I0_fit_wo = [];
    Ioffset_fit_wo = [];
    tau_fit_wo = [];
    rsq_fit_I_wo = [];
    errors = [];
    errorSTDEV = [];
    stDev_I0_wo = [];
    stDev_Ioffset_wo = [];
    stDev_tau_wo = [];
end

disp(' ') % empty line
disp('Intensity vs time exponential fit (with offset) result: ') 
disp([' I0_wo = ',num2str(I0_fit_wo),' +- ',num2str(stDev_I0_wo),';   tau_wo = ',num2str(tau_fit_wo),' +- ',num2str(stDev_tau_wo),' s.',';   Ioffset_wo = ',num2str(Ioffset_fit_wo),' +- ',num2str(stDev_Ioffset_wo),';   rsq = ',num2str(rsq_fit_I_wo)]) 

% results from exponential fit with offset of intensity vs time:
results_I_fits.I0_fit_wo = I0_fit_wo;
results_I_fits.stDev_I0_wo = stDev_I0_wo;
results_I_fits.I0_fit_percentError_wo = 100*stDev_I0_wo/I0_fit_wo;
results_I_fits.Ioffset_fit_wo = Ioffset_fit_wo;
results_I_fits.stDev_Ioffset_wo = stDev_Ioffset_wo;
results_I_fits.Ioffset_fit_percentError_wo = 100*stDev_Ioffset_wo/Ioffset_fit_wo;
results_I_fits.tau_fit_wo = tau_fit_wo;
results_I_fits.stDev_tau_wo = stDev_tau_wo;
results_I_fits.tau_fit_percentError_wo = 100*stDev_tau_wo/tau_fit_wo;
results_I_fits.rsq_fit_I_wo = rsq_fit_I_wo;
results_I_fits.Ioffset_relativeTo_I0_percent = 100*Ioffset_fit_wo/I0_fit_wo;

%% Plot results and fits:

errorbar(timeAxis,Inormalised_avg,Inormalised_error,'.k'); % plot with error bars.
hold on;
plot(fit_result_I,'b'); % plot exponential fit (no offset) as BLUE line.
plot(fit_result_I_wo,'r'); % plot exponential fit (with offset) as RED line.
legend('avgd data','exp fit-no offset','exp fit-with offset');
% legend('hide');
ylim([0 1]); 
xlim([0 n*0.04]);
xlabel('time from track start (s)'); 
ylabel('normalised intensity');
hold off;

