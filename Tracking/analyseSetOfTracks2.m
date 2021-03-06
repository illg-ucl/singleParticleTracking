function  analyseSetOfTracks2(minNumPointsInTrack,max_framesAwayFromTimeOrigin,colour_chosen,tau_fixed,min_rsq_fit_I,max_fit_error,prob_accept_fit,I_singleMolec,error_IsingleMolec)
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Analyse a set of tracks by opening the corresponding excel files in a
% number of chosen folders (for a given data set, see below) and using the info in them.
% Use all tracks with a number of points at least equal to
% minNumPointsInTrack, of the right colour ('top' or 'bottom' channel) and
% at most max_framesAwayFromTimeOrigin frames away from the time origin of
% the image sequence. Go track by track, fit it to an exponential with no
% offset and fixed tau (determined beforehand, it is an input) and determine
% (extrapolating from the fit) the initial intensity at the start of the track.
% Then only those fits with an rsquare value above the input threshold
% level min_rsq-fit_I, and with an avg_fit_error below the input threshold
%  max_fit_error, are collected into a vector.
% Then a histogram of the initial intensity values from those selected fits
% is plotted, which includes all the chosen tracks for a given data set.
% Histograms are also plotted for rsquare, percentage error of Initial
% intensity and avg_fit_error, for guidance. 
% A kernel density probability distribution for the initial intensities
% (coming from a sum of Gaussians around each initial intensity value, with
% widths given by their standard errors from the fit) is also plotted.
% From the mean and stdev of the initial intensity distribution,
% stoichiometry values can be determined.
%
% -- Inputs:
% - minNumPointsInTrack: min. no. of points a track must have in order to be
% used to form an average track. Default value is 25.
% - max_framesAwayFromTimeOrigin: max number of frames the start of the
% track can be from the time origin (when laser is switched on) for a given
% image sequence.
% - colour_chosen: this is a string. It has to be either 'top' or bottom', or 'any' for images with no top or bottom.
% Choose if you want to look at Green or Red spots (for dual label strains):
% The bottom half of the image corresponds to the Green channel -> 'bottom'.
% The top half of the image corresponds to the Red channel -> 'top'.
% Choose 'any' for images which are not divided into two channels.
% - tau_fixed: time constant in seconds of the exponential decay. It is a
% fixed/constant parameter of the exponential fits to the tracks.
% - min_rsq_fit_I: min rsq a fit must have to be accepted and its result
% stored. Default value is around 0.3 (it should be > 0).
% - max_fit_error: max average fit error (sqrt of mean sum of squares) for
% fit to be accepted. This is the max. deviation of the model from
% a given data point (on average), on the same units as the intensity axis
% of the plots. eg. for ATPase-GFP use something around 10000. Set it to
% Inf if you want no constraints, check the histogram-plot, find out the
% appropriate max value and then re-run with that as a limit.
% - prob_accept_fit: probability level to accept the fit to an exponential
% with fixed tau and no offset using a reduced chi-square statistic test.
% This should be set low, around 0.001, for example (the higher it is, the more fits will be rejected).
% If set to 0, the resulting reduced_chiSquare_limit = chi2inv(1,dof)/dof
% is infinity, and this constraint is removed.
% - I_singlemolec is the flu intensity of a single molecule of the fluorophore label.
% It is obtained from the analysis of in vitro single molecules of GFP or
% mCherry (fluorescent label). Write 1 for no calibration of stoichiometry.
% - error_IsingleMolec: error of I_singlemolec. Usually obtained as the
% standard deviation from a fitted gamma distribution to the distribution
% of in vitro intensity values for the fluorophores. It is not used for
% calculations, just to keep track of it together with the rest. Assume
% error_IsingleMolec, which is the width of the distribution is a
% systematic. Calibrate the plot of stoichiometry distrib probability
% without including it.
% 
% -- Outputs: 
% .png graphs are saved to a created folder for every fit. 
% Also graphs resulting from combining info from all tracks in data set are
% saved for histograms of the initial intensities, their percentage errors,
% rsquare of all fits and average fit error. Also graph for I0 prob
% density. These are saved in folder "results_allTracks" within the output
% folder. 
%
% The following output variables are also saved to that folder and to an excel file:
%
% Output variable "allTrack_inputs" with fields:
% data_set_label
% folder_path
% image_numbers
% minNumPointsInTrack
% max_framesAwayFromTimeOrigin
% colour_chosen
% tau_fixed
% min_rsq_fit_I
% max_fit_error
% prob_accept_fit
% num_tracks_used (final no. of accepted tracks)
%
% Output variable "results" is a structure array with fields:
%     NumDataPoints
%     FrameDistanceToTrackOrigin
%     I0_fit
%     stDev_I0
%     I0_fit_percentError
%     rsq_fit_I
%     sumOfSquares
%     sumOfSquares_normalised
%     avg_fit_error
%     tau_fixed
%     I_approx_errorBar
%     minNumPointsInTrack
%     max_framesAwayFromTimeOrigin
%     colour_chosen
%     min_rsq_fit_I
%     max_fit_error
%     prob_accept_fit
%     reduced_chiSquare
%     reduced_chiSquare_limit
%     dof
%     cumul_chi2_prob
%
% Also graphs and variables of results from all tracks together are saved in allTrack_results folder,
% and excel files of accepted tracks are re-saved in folder with extra info.  
% Output variable "allTrack_results1" is a structure where each field is a column vector. Fields are: 
%                      Initial_I
%              Initial_I_ErrorFit
%       Initial_I_percentErrorFit
%                   all_rsq_fit_I
%                  all_fit_errors
%           all_reduced_chiSquare
%     all_reduced_chiSquare_limit
%                         all_dof
%             all_cumul_chi2_prob
% 
% Output variable "allTrack_results2" is a structure with two fields as columns:
%     vector_I_initial
%      I0_prob_density
%
% The original excel files which are accepted as good fits are copied into the output folder with "bis" at the end. 
% The numbers in all fields from "results" are written into an additional
% sheet, called "initial_I data", onto those copied excel files.
%
% Example of how to run this function: 
% You need to be in the directory where the folders with the track files
% for all images in a data set are. Eg: 'C:\Isabel\ExperimData\HeikoData\ATPase-GFP\ATPase-GFP TIRF\good'.
% analyseSetOfTracks2(5,500,'bottom',2.85,0,Inf,0.001,1); Nearly no restrictions: for tracks of at least 5
% points, and they can be as far as 500 frames from the time origin and rsq >=0, and any avg_fit_error value.
% analyseSetOfTracks2(25,50,'bottom',2.85,0.2,10000,0.001,1,1); Restrictions: tracks at
% least 25 frames long and the start of the track must be at most 50 frames
% away from the time origin of the image sequence, with rsq at least 0.2 and avg_fit_error at most 10000.
%
% ADVICE: run first with no restrictions, look at the histograms and then
% run with restrictions informed from those histograms. To run with no
% restrictions, e.g.:
% analyseSetOfTracks2(5,500,'bottom',4.4,-Inf,Inf,0,1,1);
% analyseSetOfTracks2(5,500,'bottom',2.8,0,10000,0.01,379,327)

%% Build names of folders in which to find tracks to analyse:

% For building paths in which to find track excel files:
% Uncomment one of the following: 

% SINGLE LABEL STRAINS:
% ---------------------
% short_data_set_label = 'ATPase-GFP_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [89 97 99 110 112 116 122 126 128 130 134 136];
% image_numbers_1 = [89];

% short_data_set_label = 'GFP-nuoF_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [154 156 158 160 162 164 166 168 170  174];

% short_data_set_label = 'cydB-mCherry_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [228 249 253 257 261 263 265];

% short_data_set_label = 'mCherry-sdhC_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [279 281 283 285 289 291 293 295 301 303];

% short_data_set_label = 'cyoA-mCherry_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [309 313 315 317 320 322 324 326];

% short_data_set_label = 'mCherry-nuoF_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [347 349 351 353 357 359];

% short_data_set_label = 'cyoA-mCherry_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [309 313 315 317 320 322 324 326];

% DUAL LABEL STRAINS:
% -------------------
% short_data_set_label = 'cybD-mCherry-ATPase-GFp_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [498 500 509 513 515 518 522 524];
% image_numbers_1 = [498 500]; 

% short_data_set_label = 'GFP-nuoF-mCherry-sdhC_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [533 545 547 549 551 553];

% short_data_set_label = 'cydB-mCherry-GFPuv4-nuoF_'; 
% image_numbers_1 = [427 429 431 438 442 444];

% short_data_set_label = 'cydB-mCherry-GFPuv4-nuoF_'; 
% image_numbers_1 = [472 474 478 484 486];
% ---------------

% folder_names_1 = cell(1,length(image_numbers_1)); % initialise empty cell array (row).
% for i=1:length(image_numbers_1)
%     folder_name = strcat(short_data_set_label,num2str(image_numbers_1(i)));
%     folder_names_1{1,i} = folder_name;
% end

% -----------------------
% -----------------------
% or Use this when the name of the data folder is '02', contains a zero at
% start: 
% And use the right thing around line 991 too (search for "SAVE INPUTS")!!!

% short_data_set_label = 'cydB-GFP_'; % string that labels the data set and is part of the name of folders.
% % image_numbers_1 = [02];
% image_numbers_1 = {'02' '11' '13' '15' '25' '29' '37'};

short_data_set_label = 'cyoA_mCherry_'; % string that labels the data set and is part of the name of folders.
image_numbers_1 = {'3ob' '5ob' '6ob' '8ob' '10ob'};

% short_data_set_label = 'try_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = {'498' '500'};

folder_names_1 = cell(1,length(image_numbers_1)); % initialise empty cell array (row).
for i=1:length(image_numbers_1)
    folder_name = strcat(short_data_set_label,image_numbers_1{i});
    folder_names_1{1,i} = folder_name;
end
% ------------------------

all_folder_names = folder_names_1;


% ------------------
% % Use this for combining different data sets:

% short_data_set_label = 'cydB-mCherry plus mCherry-sdhC_'; % labels output folders and results
% 
% data_set_path_1 = 'C:\Isabel\ExperimData\HeikoData\cydB-mCherry\cydB-mCherry TIRF\good\cydB-mCherry_';
% image_numbers_1 = [228 249 253 257 261 263 265];
% 
% data_set_path_1bis = 'C:\Isabel\ExperimData\HeikoData\mCherry-sdhC\mCherry-sdhC TIRF\good\mCherry-sdhC_'; % string that labels the data set and is part of the name of folders
% image_numbers_1bis = [279 281 283 285 289 291 293 295 301 303];
% 
% folder_names_1 = cell(1,length(image_numbers_1)); % initialise empty cell array (row).
% for i=1:length(image_numbers_1)
%     folder_name = strcat(data_set_path_1,num2str(image_numbers_1(i)));
%     folder_names_1{1,i} = folder_name;
% end
% 
% folder_names_1bis = cell(1,length(image_numbers_1bis)); % initialise empty cell array (row).
% for i=1:length(image_numbers_1bis)
%     folder_name = strcat(data_set_path_1bis,num2str(image_numbers_1bis(i)));
%     folder_names_1bis{1,i} = folder_name;
% end
% 
% all_folder_names = [folder_names_1 folder_names_1bis];
% --------------------------


%% Choosing channel (top or bottom, red or green channel):

if strcmp(colour_chosen,'any')
    any_colour =1;
else
    any_colour = 0;
end
    
%% Create new directory for saving plots and results:

% Make new folder (new directory) to save trajectory analysis results:
initial_folder_path = cd;
% Create new folder for outputs inside current folder:
output_folder_name = strcat(short_data_set_label,'I0_',colour_chosen);
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(output_folder_name); % make new directory.

output_folder_path = strcat(initial_folder_path,'\',output_folder_name);
% Create another folder inside directory output_folder_path, in which to
% save results from combining all tracks:
output_folder_name2 = 'results_allTracks';
cd(output_folder_path);
mkdir(output_folder_name2); % make new directory.
cd(initial_folder_path);


%% Go through track files (xls files) in each of the previous folders:

% Initialise index of accepted tracks to accummulate for further analysis:
m = 1;
% This is the number of tracks with at least minNumPointsInTrack points in
% them, which are close enough to the time origin and which are 'top' or
% 'bottom' as requested through input.

for j =1:length(all_folder_names) % j is the index for the number of folders to go through.
    
    cd(all_folder_names{j}); % move into folder (directory);
    
    xlsFileNames0 = dir('*.xls'); % Get names of excel files in that folder (track analysis .xls files).
    % Error control:
    if isempty(xlsFileNames0) % If there are no .xls track data files, show error and exit function:
        error('Check you are in the correct directory. No .xls files found in folder.');
    end
    
    xlsFileNames = {xlsFileNames0.name}; % cell array of strings.
    
    % Loop through each track-analysis xls-file inside folder j:
    for k=1:length(xlsFileNames) % k is the index for each excel file.
        
        disp(xlsFileNames{k});
        
        pos1 = strfind(xlsFileNames{k},'.xls'); % position of the start of the string '.xls' in the excel file name.
        figName_0 = xlsFileNames{k}(1:(pos1-1)); % Take the name of the input excel file (with the end bit '.xls' removed) as the new name for saving a later .png figure.
        figName = strcat(figName_0,'_FitFor_Initial_I'); % name of figure file to save to.
        
        % Open and read the data:
        % -----------------------
        % Import the data in the sheet named 'Track results':
        [numeric,txt,raw] = xlsread(xlsFileNames{k},'Track info');
        % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
        str_TrackInfo = cell2struct(raw(:,2),raw(:,1),1);
        NumDataPoints = str_TrackInfo.NumDataPoints;
        FrameForTimeOrigin = str_TrackInfo.FrameForTimeOrigin; % frame for time origin for a given data set (frame for which laser is switched on).
        FirstTrajFrame = str_TrackInfo.FirstTrajFrame; % First frame in this particular track.
        AbsTimeOrigin = str_TrackInfo.AbsTimeOrigin; % time origin for a given data set (frame for which laser is switched on).
        % TrajStartTime = str_TrackInfo.TrajStartTime; % time at the start of a particular track.
        TimeBetweenFrames = str_TrackInfo.TimeBetweenFrames; % in seconds.
        TopOrBottom = str_TrackInfo.TopOrBottom; % 'top' or 'bottom' region of image (colour channel, red or green);
        FrameDistanceToTrackOrigin = FirstTrajFrame-FrameForTimeOrigin;        
        % -----------------------
        % Import the data in the sheet named 'Mobility results':
        [numeric2,txt2,raw2] = xlsread(xlsFileNames{k},'Mobility results');
        % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
        str_Mobility = cell2struct(raw2(:,2),raw2(:,1),1);
        % Later we save the following fields from strMobility as fields of
        % results(m):
        %         fit_msd_line_rsq
        %         fit_msd_conf_rsq
        %         diffusion1D_coeff_micronsSqrdPerSec
        %         diffusion1D_coeff_error
        %         size_1Dconfin_nm
        %         size_1Dconfin_nm_error
        %         Brownian_flag
        %         Confined_flag
        %         OtherMobility_flag
        % -----------------------
        % Import the data in the sheet named 'cell coordinates':
        [numeric3,txt3,raw3] = xlsread(xlsFileNames{k},'cell coordinates');
        % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
        str_cellCoords = cell2struct(raw3(:,2),raw3(:,1),1);
        % Later we save the following fields from str_cellCoords as fields of
        % results(m):
        %         cell_width_nm
        %         cell_length_nm
        %         traj_in_pole

        
        % -----
        % Use a criterium to accept which track data we accummulate to use it later:
        if NumDataPoints >= minNumPointsInTrack && ... % if num of points in track is at least minNumPointsInTrack (25):
                (any_colour ==1 || strcmp(TopOrBottom,colour_chosen) == 1) && ...% if the track is of the right colour/channel
                FrameDistanceToTrackOrigin <= max_framesAwayFromTimeOrigin % if the track is close enough to the time origin of the image sequence.            
            
            % Import the data in the sheet named 'Results I fits':
            [numeric,txt,raw] = xlsread(xlsFileNames{k},'Results I fits');
            % Turn imported data from excel file into a structure:
            str_ResultsIfits = cell2struct(raw(:,2),raw(:,1),1);
            I0_fit_prelim = str_ResultsIfits.I0_fit; % preliminary value used for guessing, from fit with no offset.
            % Now the fit with no offset (used later to find an approx
            % error bar for data points):
            I0_fit_wo_prelim = str_ResultsIfits.I0_fit_wo;
            Ioffset_fit_wo_prelim = str_ResultsIfits.Ioffset_fit_wo;
            tau_fit_wo_prelim = str_ResultsIfits.tau_fit_wo;                     
            
            % Import the data in the sheet named 'Track data':
            [numeric,txt,raw] = xlsread(xlsFileNames{k},'Track data');
            % Turn imported data from excel file into a structure:
            str_TrackData = struct; % create empty structure to fill up.
            for i = 1:length(txt)
                str_TrackData = setfield(str_TrackData, txt{i}, numeric(:,i)); % create field in the structure.
                % Each field in the structure contains a column vector with the data.
            end
            
            % Vectors useful for later:
            % Time in seconds relative to the time origin of the image sequence:
            time_rel_to_origin = str_TrackData.timeabs-AbsTimeOrigin; % column vector
            % integrated spot intensity:
            intensity = str_TrackData.intensity; % column vector            
            
            % Flatten intensity data by subtracting its fit to an
            % exponential with offset, to then get an approx. average measure of
            % the error bar of the original track data from the standard
            % deviation of the flattened data, remaining I_flattened:
            I_flattened = intensity - (Ioffset_fit_wo_prelim + I0_fit_wo_prelim*exp(-time_rel_to_origin/tau_fit_wo_prelim));
            % consider the standard deviation of the remaining noise as the
            % approx. error bar of the original track data:
            I_approx_errorBar = std(I_flattened); 
            % % Check graphically:
            % plot(time_rel_to_origin,I_flattened)
            
            % ---------------------------------------------------
            % Fit intensity data to an exponentially decaying function (with no offset):
            
            exp_no_offset = fittype('I0*exp(-t/tau)','independent','t','problem','tau'); % define exponential funtion to fit to, with 't' as independent
            % variable and 'tau' as a fixed parameter (constant).
            options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
            % Guesses for fit parameters:
            guess_I0 = I0_fit_prelim; % since I has been normalised.
            % Use coeffnames(fit_result_I) later to find out order of parameters.
            options.StartPoint = guess_I0; % give guess parameters for fit. This avoids a warning message. Give in right order!.
            options.Lower = 0; % Lower bounds for fit parameters. In order: I0.
            options.Upper = 10*I0_fit_prelim; % Upper bounds for fit parameters. In order: I0.
            [fit_result_I gof] = fit(time_rel_to_origin,intensity,exp_no_offset,options,'problem',tau_fixed); % fit_result_I contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
            % fit_param_names = coeffnames(fit_result_I); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
            fit_param_values = coeffvalues(fit_result_I); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
            I0_fit = fit_param_values(1); % I0 intensity value from fit.          
            rsq_fit_I = gof.rsquare; % rsquare coefficient of fit.
            sumOfSquares = gof.sse;
            sumOfSquares_normalised = sumOfSquares/length(intensity); % normalised dividing by the no. of points in trace to get the average squared deviation of model from data.
            avg_fit_error = sqrt(sumOfSquares_normalised); % approx average deviation of a data point from the fitted model.
            errors = confint(fit_result_I,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
            errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
            stDev_I0 = errorSTDEV(1);
            % Error control, if the error is NaN (not a number):
            if isnan(stDev_I0)
                % Approximate the error bar of I0 with the approximate
                % experimental error bar:
                stDev_I0 = I_approx_errorBar;
            end
            
            % Calculate reduced Chi-square test goodness of fit:
            % Consider the error bar is the same for all points in track
            % and equal to I_approx_errorBar and calculate the chi square:
            chiSquare_fit = sumOfSquares/(I_approx_errorBar^2);
            nparams_fit = 1; % number of free parameters in the fit, just I0.
            dof = length(intensity)-nparams_fit-1; % number of degrees of freedom.
            reduced_chiSquare = chiSquare_fit/dof; % reduced chi squared statistic. Should be equal to 1 for a good fit.
            % Maximum value the reduced_chiSquare can have for the fit to be accepted,
            % with a rather low probability level of prob_accept_fit 0.001 (0.1%):
            reduced_chiSquare_limit = chi2inv(1-prob_accept_fit,dof)/dof; % inverse chi-square cumulative distribution function. 
            % All values in the chi-square distrib with dof degrees of freedom will be below reduced_chiSquare_limit with a probability of 1-prob_accept_fit. 
            % If prob_accept_fit = 0, reduced_chiSquare_limit = Inf, so the
            % constraint is removed.
            % Prob. that a value in the chi-square distrib with dof degrees
            % of freedom is below chiSquare_fit:
            cumul_chi2_prob = chi2cdf(chiSquare_fit,dof);
            % [reduced_chiSquare, reduced_chiSquare_limit dof cumul_chi2_prob]
            
            % Store result only if it fulfills certain conditions:
            % These are the final accepted tracks.
            if rsq_fit_I >= min_rsq_fit_I && ...
                    avg_fit_error <= max_fit_error && ...
                    reduced_chiSquare <= reduced_chiSquare_limit
                
                %             disp(' ') % empty line
                %             disp('Intensity vs time exponential fit (with no offset) and tau fixed: ')
                %             disp([' I0 = ',num2str(I0_fit),' +- ',num2str(stDev_I0),';  tau_fixed = ',num2str(tau_fixed),';   rsq = ',num2str(rsq_fit_I)])
                
                % Stoichiometry determination: number of molecules:
                stoichiom = round(I0_fit/I_singleMolec); % I_singlemolec is the flu intensity of a single molecule of the fluorophore label.
                
                % Save to output: results from exponential fit with no offset and fixed tau, of intensity vs time:
                results(m).NumDataPoints = NumDataPoints;
                results(m).FrameDistanceToTrackOrigin = FrameDistanceToTrackOrigin;
                results(m).Iavg = mean(intensity); % average intensity of data set (used for short tracks later on).
                results(m).I0_fit = I0_fit;
                results(m).stDev_I0 = stDev_I0;
                results(m).I0_fit_percentError = 100*stDev_I0/I0_fit;
                results(m).rsq_fit_I = rsq_fit_I;
                results(m).I_singleMolec = I_singleMolec; % save input to results. I_singlemolec is the flu intensity of a single molecule of the fluorophore label.           
                results(m).error_IsingleMolec = error_IsingleMolec;  % save input to results. Error (std) of the single molec intensity value, from std of the gamma distrib fitted to histogram of intensities of in vitro data. 
                results(m).stoichiom = stoichiom; % stoichiometry, number of molecules. Rounded to nearest integer.     
                results(m).error_stoichiom = round(stDev_I0/I_singleMolec); % error for plot. Assume error_IsingleMolec, which is the width of the distribution is a systematic. Calibrate the plot without including it.
                % results(m).error_stoichiom = round((I0_fit/I_singleMolec)*sqrt((stDev_I0/I0_fit)^2+(error_IsingleMolec/I_singleMolec)^2));
                results(m).sumOfSquares = sumOfSquares;
                results(m).sumOfSquares_normalised = sumOfSquares_normalised; % divided by no. of points in track
                results(m).avg_fit_error = avg_fit_error;
                results(m).tau_fixed = tau_fixed; % save input
                results(m).I_approx_errorBar = I_approx_errorBar;
                results(m).minNumPointsInTrack = minNumPointsInTrack; % save input
                results(m).max_framesAwayFromTimeOrigin = max_framesAwayFromTimeOrigin; % save input
                results(m).colour_chosen = colour_chosen; % save input
                results(m).min_rsq_fit_I = min_rsq_fit_I; % save input 
                results(m).max_fit_error = max_fit_error; % save input 
                results(m).prob_accept_fit = prob_accept_fit; % save input 
                results(m).reduced_chiSquare = reduced_chiSquare;
                results(m).reduced_chiSquare_limit = reduced_chiSquare_limit;
                results(m).dof = dof;
                results(m).cumul_chi2_prob = cumul_chi2_prob;
                
                % Save directly to results useful results from track excel
                % file, from several "Sheets":
                % From sheet "Track info":
                results(m).Track_meanXvalue_cellCoords = str_TrackInfo.Track_meanXvalue_cellCoords;
                results(m).Track_meanYvalue_cellCoords = str_TrackInfo.Track_meanYvalue_cellCoords;
                results(m).meanSpotWidth = str_TrackInfo.meanSpotWidth;
                % From sheet "Mobility results":
                results(m).fit_msd_line_rsq = str_Mobility.fit_msd_line_rsq;
                results(m).fit_msd_conf_rsq = str_Mobility.fit_msd_conf_rsq;
                results(m).diffusion1D_coeff_micronsSqrdPerSec = str_Mobility.diffusion1D_coeff_micronsSqrdPerSec;
                results(m).diffusion1D_coeff_error = str_Mobility.diffusion1D_coeff_error;
                results(m).size_1Dconfin_nm = str_Mobility.size_1Dconfin_nm;
                results(m).size_1Dconfin_nm_error = str_Mobility.size_1Dconfin_nm_error;
                results(m).Brownian_flag = str_Mobility.Brownian_flag;
                results(m).Confined_flag = str_Mobility.Confined_flag;
                results(m).OtherMobility_flag = str_Mobility.OtherMobility_flag;
                % From sheet "cell coordinates":
                results(m).cell_width_nm = str_cellCoords.cell_width_nm;
                results(m).cell_length_nm = str_cellCoords.cell_length_nm;
                results(m).traj_in_pole = str_cellCoords.traj_in_pole;            

                        
                % ---------------------------------------------------
                % Plot results and fits:
                % Create figure. 'position' vector is [left, bottom, width, height].
                h1 = figure('units','inches','position',[8 4 6 5]);
                plot(time_rel_to_origin,intensity,'.k'); % plot intensity trace.
                hold on;
                full_time_axis = (0:TimeBetweenFrames:500*TimeBetweenFrames);
                plot(full_time_axis,I0_fit*exp(-full_time_axis/tau_fixed),'-b'); % plot exponential fit (no offset) as BLUE line.
                legend('data','exp fit-no offset');
                % legend('hide');
                % ylim([0 1.2]);
                % xlim([0 500*TimeBetweenFrames]); % usually there are no more than 500 frames
                xlabel('time relative to origin (s)');
                ylabel('integrated spot intensity');
                hold off;
                
                str1(1) = {[xlsFileNames{k}]};
                str1(2) = {['No. data points in trajectory: ',num2str(NumDataPoints)]};
                str1(3) = {['First frame in trajectory:  ',num2str(FirstTrajFrame)]};
                str1(4) = {['Frame for sequence time origin: ',num2str(FrameForTimeOrigin)]};
                str1(5) = {['Time between frames (seconds): ',num2str(TimeBetweenFrames)]};                
                str1(6) = {['colour_chosen: ',num2str(colour_chosen)]};
                str1(7) = {['minNumPointsInTrack: ',num2str(minNumPointsInTrack)]};
                str1(8) = {['max_framesAwayFromTimeOrigin: ',num2str(max_framesAwayFromTimeOrigin)]}; 
                str1(9) = {['min_rsq_fit_I: ',num2str(min_rsq_fit_I)]}; 
                str1(10) = {['max_fit_error: ',num2str(max_fit_error)]};                 
                str1(11) = {[' I0_expFit = ',num2str(I0_fit),' +- ',num2str(stDev_I0)]};
                str1(12) = {[' tau_fixed = ',num2str(tau_fixed),' s']};
                str1(13) = {[' rsq_fit = ',num2str(rsq_fit_I)]};
                str1(14) = {[' avg fit error = ',num2str(avg_fit_error)]};
                str1(15) = {[' cumul_chi2_prob = ',num2str(cumul_chi2_prob)]};
                str1(16) = {[' stoichiom = ',num2str(stoichiom)]};
                % Display text (text(x,y,'string','PropertyName',PropertyValue....)):
                axes_limits = axis; % get axes limits of current figure: axes_limits = [xmin xmax ymin ymax].
                % ylim([0 axes_limits(4)]); % make the lower y-axis limit be zero.
                text(axes_limits(2)/3,axes_limits(4)/1.8,str1,'Interpreter','none')  % display path of image sequence at top left of figure at position (xmin,ymax) of axes.
                % doing: ...'Interpreter','latex' creates subindices for _staff...
                
                % ---------------------------------------------------     
                % Copy the original excel file for the ACCEPTED track
                % (xlsFileNames{k}) into the output folder,
                % copyfile('source','destination'): 
                excelfilename_0 = xlsFileNames{k}; % original excel file name, including .xls.
                excelfilename_1 = excelfilename_0(1:findstr(excelfilename_0,'.xls')-1); % original excel file name, without the file extension .xls.
                excelfilename_2 = strcat(excelfilename_1,'_bis','.xls'); % Add '_bis' to distinguish it from original excel file.
                track_excel_file_newpath = strcat(output_folder_path,'\',excelfilename_2);
                copyfile(xlsFileNames{k},track_excel_file_newpath)
                     
                % Update the track excel file adding new results from the fit for I0 to a new sheet (including stoichiometry value):
                dataForSheet_I0 = [fieldnames(results(m)) struct2cell(results(m))]; % numbers, first column is fieldnames of element {m}, second one is values.
                warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.
                xlswrite(track_excel_file_newpath,dataForSheet_I0,'initial_I data'); % write data to sheet 'initial_I data' in excel file.
                % ----------------------------------------------------
                
                % cd('..'); % go back to previous directory where folder has been created to save images.
                cd(initial_folder_path) % go to initial directory where folder has been created to save output.
                % Save figure as png in folder created (output_folder_name), and close figure:
                saveFigurePNG(output_folder_name,figName)
                % move back into folder (directory) where excel files are:
                cd(all_folder_names{j});
                
                m = m+1; % advance index for accummulating results.
                
            end  % if rsq_fit_I >= min_rsq_fit_I
            
        end % if NumDataPoints >= minNumPointsInTrack.
        
    end
    
    % cd('..'); % go back to previous directory.
     cd(initial_folder_path) % go to initial directory where folder has been created to save output.
     
end

num_tracks_used = m-1;

disp(' '); % empty line.
disp(['The number of tracks fulfilling conditions is: ',num2str(num_tracks_used),'']);

% -------------------------------
% All-track results, column vectors (for output):
allTrack_results1.NumDataPoints = [results.NumDataPoints]'; % no. of data points in track.
allTrack_results1.FrameDistanceToTrackOrigin = [results.FrameDistanceToTrackOrigin]'; 
allTrack_results1.Initial_I = [results.I0_fit]'; 
allTrack_results1.Initial_I_ErrorFit =[results.stDev_I0]'; % error from fit, standar deviation.
allTrack_results1.Initial_I_percentErrorFit = [results.I0_fit_percentError]'; % percentage error from fit
allTrack_results1.all_stoichiom = [results.stoichiom]'; 
allTrack_results1.all_error_stoichiom = [results.error_stoichiom]'; 
allTrack_results1.all_rsq_fit_I = [results.rsq_fit_I]';
allTrack_results1.all_fit_errors = [results.avg_fit_error]'; % avg deviation of fitted model from data.
allTrack_results1.all_reduced_chiSquare = [results.reduced_chiSquare]';
allTrack_results1.all_reduced_chiSquare_limit =[results.reduced_chiSquare_limit]';
allTrack_results1.all_dof = [results.dof]';
allTrack_results1.all_cumul_chi2_prob = [results.cumul_chi2_prob]';     
allTrack_results1.all_Iavg = [results.Iavg]'; % Mean intensity of all points in track (used later on for short tracks).
allTrack_results1.all_I_approx_errorBar = [results.I_approx_errorBar]'; % std of flattened I trace (subtracting fit to exp with offset from data).
allTrack_results1.all_I_approx_errorBar_sqr = (allTrack_results1.all_I_approx_errorBar).^2; % same as previous, squared.
          
allTrack_results1.all_Track_meanXvalue_cellCoords = [results.Track_meanXvalue_cellCoords]';
allTrack_results1.all_Track_meanYvalue_cellCoords = [results.Track_meanYvalue_cellCoords]';
allTrack_results1.all_meanSpotWidth = [results.meanSpotWidth]';
allTrack_results1.all_fit_msd_line_rsq = [results.fit_msd_line_rsq]';
allTrack_results1.all_fit_msd_conf_rsq = [results.fit_msd_conf_rsq]';
allTrack_results1.all_diffusion1D_coeff_micronsSqrdPerSec = [results.diffusion1D_coeff_micronsSqrdPerSec]';
allTrack_results1.all_diffusion1D_coeff_error = [results.diffusion1D_coeff_error]';
allTrack_results1.all_size_1Dconfin_nm = [results.size_1Dconfin_nm]';
allTrack_results1.all_size_1Dconfin_nm_error = [results.size_1Dconfin_nm_error]';
allTrack_results1.all_Brownian_flag = [results.Brownian_flag]';
allTrack_results1.all_Confined_flag = [results.Confined_flag]';
allTrack_results1.all_OtherMobility_flag = [results.OtherMobility_flag]';
allTrack_results1.all_cell_width_nm = [results.cell_width_nm]';
allTrack_results1.all_cell_length_nm = [results.cell_length_nm]';
allTrack_results1.all_traj_in_pole = [results.traj_in_pole]';

%---------------------------------
% Creatte plot with histograms no.2:
% Histograms of initial intensity for all tracks (with accepted fits) and of its error, and of
% stoichiometry and of its error.
figure('units','inches','position',[6 3 8 8]);

subplot(2,2,1)
hist(allTrack_results1.Initial_I,50) % 50 bins in histogram
xlabel('Initial Intensity');
ylabel('frequency');
axes_limits_I = axis; % get axes limits of current figure as [xmin xmax ymin ymax].
% Get limits of initial intensity axis of histogram:
I_initial_min = 0;
I_initial_max = axes_limits_I(2);

subplot(2,2,2)
hist(allTrack_results1.Initial_I_percentErrorFit,50) % 50 bins in histogram
xlabel('Initial-I percentError');
ylabel('frequency');

subplot(2,2,3)
hist(allTrack_results1.all_stoichiom,50) % 50 bins in histogram
xlabel('stoichiometry (no. molecs)');
ylabel('frequency');
axes_limits_stoich = axis; % get axes limits of current figure as [xmin xmax ymin ymax].
% Get limits of stoichiom axis of histogram:
stoichiom_min = 0;
stoichiom_max = axes_limits_stoich(2);

subplot(2,2,4)
hist(allTrack_results1.all_error_stoichiom,50) % 50 bins in histogram
xlabel('error-stoichiom');
ylabel('frequency');

% Save current figure of histograms from all tracks as .png:
cd(output_folder_path);
figName1 = strcat(short_data_set_label,colour_chosen,'_Histograms'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName1)
cd(initial_folder_path);

%---------------------------------
% Creatte plot with histograms no.2:
% Histograms of errors and goodness of fit params used to accept or reject
% fits. Histograms for all tracks (with accepted fits) of rsq_fit,
% avg_fit_error and reduced_chiSquare.
figure('units','inches','position',[6 3 8 8]);

subplot(2,2,1)
hist(allTrack_results1.all_rsq_fit_I,50) % 50 bins in histogram
xlabel('rsquare of fits (tau fixed)');
ylabel('frequency');

subplot(2,2,2)
hist(allTrack_results1.all_fit_errors,50) % 50 bins in histogram
xlabel('avg fit error');
ylabel('frequency');

subplot(2,2,3)
hist(allTrack_results1.all_reduced_chiSquare,50) % 50 bins in histogram
xlabel('reduced chiSquare');
ylabel('frequency');

subplot(2,2,4)
hist(allTrack_results1.all_cumul_chi2_prob,50) % 50 bins in histogram
xlabel('cumul chi2 prob');
ylabel('frequency');

% Save current figure of histograms from all tracks as .png:
cd(output_folder_path);
figName2 = strcat(short_data_set_label,colour_chosen,'_FitTestHistogr'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName2)
cd(initial_folder_path);


%---------------------------------
% Plot PROBABILITY DENSITY of I0 values (initial intensities):
% Create plot with kernel density plot (sum of Gaussians to generate a
% probability density instead of a histogram):
% Kernel density probability distribution for the initial intensities
% (coming from a sum of Gaussians around each initial intensity value, with
% widths given by their standard errors from the fit).

% Create normalised Gaussians (one around each initial intensity value) and
% add them up:
% Vector covering the whole range of initial intensities with 10^6
% points (x vector for gaussians). Careful! Need high enough resolution for this, so large no. of points!:
vector_I_initial = (I_initial_min:I_initial_max/1000000:I_initial_max)'; % column vector
% Initialise empty vector in which to accumulate normalised Gaussians:
vector_all_Gaussians = zeros(length(vector_I_initial),1); % column vector

% Loop through all accepted fits of tracks:
for i = 1:length(allTrack_results1.Initial_I)
    
    % normalised Gaussian centered around I_initial value and error indexed
    % by i:
    normalised_gaussian_vector = normalisedGaussian(vector_I_initial,allTrack_results1.Initial_I(i),allTrack_results1.Initial_I_ErrorFit(i));
    % Error control in case any element is NaN:
    if ~isempty(normalised_gaussian_vector(isnan(normalised_gaussian_vector))) % if there are any NaN elements
        for j = 1:length(normalised_gaussian_vector) % convert all NaN elements into zeros:
            if isnan(normalised_gaussian_vector(j))
                normalised_gaussian_vector(j) = 0;
            end
        end

    end
    % plot(vector_I_initial,normalised_gaussian_vector)
    
    % Accummulate all normalised Gaussians onto final vector:
    vector_all_Gaussians = vector_all_Gaussians + normalised_gaussian_vector;
end

% Normalise vector which accummulates all Gaussians by dividing by the
% number of Gaussians (number of accepted track fits), and then plot:
% The total area under the curve should be 1.
I0_prob_density = vector_all_Gaussians/length(allTrack_results1.Initial_I);
plot(vector_I_initial,I0_prob_density)
xlabel('Initial Intensity');
ylabel('Probability density');
% Save current figure and close it (I0 prob density from all tracks) as .png in 'results_allTracks' folder:
cd(output_folder_path);
figName3 = strcat(short_data_set_label,colour_chosen,'_InitialI_ProbDensity0'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName3)

% Plot again and hold:
plot(vector_I_initial,I0_prob_density)
xlabel('Initial Intensity');
ylabel('Probability density');
hold on; % hold it to plot the fit to a gamma distrib afterwards.

% Resample prob density vs I0 data by interpolating, to reduced the no. of
% points for saving:
x1 = vector_I_initial; % column vector of original x data.
y1 = I0_prob_density;  % column vector of original y data.
x1min = min(x1);
x1max = max(x1);
x2 = (x1min:((x1max-x1min)/999):x1max)'; % new x vector with 1000 points, for re-sampling.
y2 = interp1q(x1,y1,x2); % interpolated prob density: one point for each value of x2.
% Resampled vectors for I0 prob density plot:
vector_I_initial_resampled = x2;
I0_prob_density_resampled = y2; 

% All-track results, save the probability density data as output:
allTrack_results2.I0_axis = vector_I_initial_resampled; % column vector
allTrack_results2.I0_prob_density = I0_prob_density_resampled; % column vector


% ----------------
% Fit histogram prob density of I0 values to a Gamma distribution:

% Preliminarly find the parameters of the gamma distribution which best
% fits the vector of all I0 values and use later as guesses for
% fitting the prob density distrib of I0 values:
[params_gamma, pci] = gamfit(allTrack_results1.Initial_I);
% For a gamma distribution with parameters a and b, gampdf(x,a,b), params_gamma(1)
% gives "a" and params_gamma(2) gives "b". "pci" gives the 95% confidence intervals
% for "a" and "b" in each column, respectively. The 95% confidence intervals for "a" is the first column, pci(1,1) to
% pci(2,1), and for "b" it is the second column: pci(1,2) to pci(2,2).

fun_to_fit = fittype('gampdf(x,a,b)*c','independent','x'); % define exponential funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
options.StartPoint = [params_gamma(1) params_gamma(2) 1]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
% options.Lower = [0 0.1*params_gamma(2) 0]; % give lower bound for parameters. 
% options.Upper = [10 10*params_gamma(2) 2]; % give upper bound for parameters. 
% c: guess is 1, since it should already be normalised.
try
    [fit_result gof] = fit(vector_I_initial_resampled,I0_prob_density_resampled,fun_to_fit,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
    % fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
    a_fit = fit_param_values(1);
    b_fit = fit_param_values(2);
    c_fit = fit_param_values(3);
    rsq_fit = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    stDev_a = errorSTDEV(1);
    stDev_b = errorSTDEV(2);
    stDev_c = errorSTDEV(3);
      
    % disp(' ') % empty line
    % disp('Fit of I0 prob density to gamma distrib: ')
    % disp([' a = ',num2str(a_fit),' +- ',num2str(stDev_a),';   b = ',num2str(b_fit),' +- ',num2str(stDev_b)])
    % disp([' c = ',num2str(c_fit),' +- ',num2str(stDev_c),';    rsqr = ',num2str(rsq_fit)])
    
    % Add plot of fit to previously held plot:
    plot(fit_result,'-r'); % plot fit to a gamma distrib as a red curve
    xlabel('I0');
    ylabel('Probability density')
    
catch ME1 % if the fit failed (shape can be not gamma distrib-like)
    a_fit = params_gamma(1);
    b_fit = params_gamma(2);
    c_fit = 0;
    rsq_fit = 0;
    stDev_a = (pci(2,1)-pci(1,1))/2;
    stDev_b = (pci(2,2)-pci(1,2))/2;
    stDev_c = 0;
end
    
hold off; % for figure
    
[mean_gamma,variance_gamma] = gamstat(a_fit,b_fit); % returns mean of and variance for the gamma distribution with parameters specified by params_gamma(1) and params_gamma(2).
% disp('Gamma stats from new fit: ')
% disp(['params gamma distrib:  ',num2str([a_fit b_fit])]);
% disp(['mean_gamma: ',num2str(mean_gamma)]);
% disp(['std_gamma: ',num2str(sqrt(variance_gamma))]);
% disp(['mode gamma: ',num2str((a_fit-1)*b_fit)]);

% numbers for output:
allTrack_I0.a_fit = a_fit;
allTrack_I0.error_a_fit = stDev_a;
allTrack_I0.b_fit = b_fit;
allTrack_I0.error_b_fit = stDev_b;
allTrack_I0.c_fit = c_fit;
allTrack_I0.error_c_fit = stDev_c;
allTrack_I0.rsq_fit = rsq_fit;
allTrack_I0.mean_gamma = round(mean_gamma);
allTrack_I0.mode_gamma_fit = round((a_fit-1)*b_fit); % I0 value (x-value) corresponding to the maximum of the fitted curve.
allTrack_I0.std_gamma = round(sqrt(variance_gamma));
allTrack_I0.mode_data = round(vector_I_initial(I0_prob_density==max(I0_prob_density))); % I0 value (x-value) corresponding to the maximum of the data curve.

% ---

% Save current figure (I0 prob density from all tracks) as .png in 'results_allTracks' folder:
cd(output_folder_path);
figName3 = strcat(short_data_set_label,colour_chosen,'_InitialI_ProbDensity'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName3)


%---------------------------------
%---------------------------------
% STOICHIOMETRY:
% Plot PROBABILITY DENSITY of STOICHIOMETRY values:
% Create plot with kernel density plot (sum of Gaussians to generate a
% probability density instead of a histogram):
% Kernel density probability distribution for the stoichiometry values
% (sum of Gaussians around each stoichiom value, with widths given by their errors).

% Create normalised Gaussians (around each stoichiom value) and add them up:
% Vector covering the whole range of initial intensities with 10^6
% points (x-vector for gaussians). Careful! Need high enough resolution for this, so large no. of points!:
vector_stoichiom = (stoichiom_min:stoichiom_max/1000000:stoichiom_max)'; % column vector
% Initialise empty y-vector in which to accumulate normalised Gaussians:
vector_all_Gaussians_stoich = zeros(length(vector_stoichiom),1); % column vector

% Loop through all accepted fits of tracks:
for i = 1:length(allTrack_results1.all_stoichiom)
    
    % normalised Gaussian centered around stoichiom value and error indexed by i:
    normalised_gaussian_vector_stoich = normalisedGaussian(vector_stoichiom,allTrack_results1.all_stoichiom(i),allTrack_results1.all_error_stoichiom(i));
    % Error control in case any element is NaN:
    if ~isempty(normalised_gaussian_vector_stoich(isnan(normalised_gaussian_vector_stoich))) % if there are any NaN elements
        for j = 1:length(normalised_gaussian_vector_stoich) % convert all NaN elements into zeros:
            if isnan(normalised_gaussian_vector_stoich(j))
                normalised_gaussian_vector_stoich(j) = 0;
            end
        end

    end
    % plot(vector_I_initial,normalised_gaussian_vector)
    
    % Accummulate all normalised Gaussians onto final vector:
    vector_all_Gaussians_stoich = vector_all_Gaussians_stoich + normalised_gaussian_vector_stoich;
end

% Normalise vector which accummulates all Gaussians by dividing by the
% number of Gaussians (number of accepted track fits), and then plot:
% The total area under the curve should be 1.
stoich_prob_density = vector_all_Gaussians_stoich/length(allTrack_results1.all_stoichiom);
% Plot:
plot(vector_stoichiom,stoich_prob_density) % curve in blue
xlabel('no. of molecules (stoichiometry)');
ylabel('Probability density');
% Save current figure and close it (stoichiometry prob density from all tracks) as .png in 'results_allTracks' folder:
cd(output_folder_path);
figName4 = strcat(short_data_set_label,colour_chosen,'_stoichiom_ProbDensity0'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName4)

% Plot again and hold:
plot(vector_stoichiom,stoich_prob_density) % curve in blue
xlabel('no. of molecules (stoichiometry)');
ylabel('Probability density');
hold on; % hold it to plot the fit to a gamma distrib afterwards.

% Resample prob density vs stoichiometry data by interpolating, to reduced the no. of
% points for saving:
x1 = vector_stoichiom; % column vector of original x data.
y1 = stoich_prob_density;  % column vector of original y data.
x1min = min(x1);
x1max = max(x1);
x2 = (x1min:((x1max-x1min)/999):x1max)'; % new x vector with 1000 points, for re-sampling.
y2 = interp1q(x1,y1,x2); % interpolated prob density: one point for each value of x2.
% Resampled vectors for stoichiometry prob density plot:
vector_stoichiom_resampled = x2;
stoich_prob_density_resampled = y2; 

% All-track results, save the probability density data as output:
allTrack_results3.stoichiom_axis = vector_stoichiom_resampled; % column vector
allTrack_results3.stoichiom_prob_density = stoich_prob_density_resampled; % column vector


% ----------------
% Fit histogram prob density of stoichiometry values to a Gamma distribution:

% Preliminarly find the parameters of the gamma distribution which best
% fits the vector of all stoichiometry values and use later as guesses for
% fitting the prob density distrib of stoichiometry values:
[params_gamma, pci] = gamfit(allTrack_results1.all_stoichiom);
% For a gamma distribution with parameters a and b, gampdf(x,a,b), params_gamma(1)
% gives "a" and params_gamma(2) gives "b". "pci" gives the 95% confidence intervals
% for "a" and "b" in each column, respectively. The 95% confidence intervals for "a" is the first column, pci(1,1) to
% pci(2,1), and for "b" it is the second column: pci(1,2) to pci(2,2).

fun_to_fit = fittype('gampdf(x,a,b)*c','independent','x'); % define exponential funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
options.StartPoint = [params_gamma(1) params_gamma(2) 1]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
% options.Lower = [0 0.1*params_gamma(2) 0]; % give lower bound for parameters. 
% options.Upper = [10 10*params_gamma(2) 2]; % give upper bound for parameters. 
% c: guess is 1, since it should already be normalised.
try
    [fit_result gof] = fit(vector_stoichiom_resampled,stoich_prob_density_resampled,fun_to_fit,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
    % fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
    a_fit = fit_param_values(1);
    b_fit = fit_param_values(2);
    c_fit = fit_param_values(3);
    rsq_fit = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    stDev_a = errorSTDEV(1);
    stDev_b = errorSTDEV(2);
    stDev_c = errorSTDEV(3);
      
    % disp(' ') % empty line
    % disp('Fit of stoichiom prob density to gamma distrib: ')
    % disp([' a = ',num2str(a_fit),' +- ',num2str(stDev_a),';   b = ',num2str(b_fit),' +- ',num2str(stDev_b)])
    % disp([' c = ',num2str(c_fit),' +- ',num2str(stDev_c),';    rsqr = ',num2str(rsq_fit)])
    
    % Add plot of fit to previously held plot:
    plot(fit_result,'-r'); % plot fit to a gamma distrib as a red curve
    xlabel('no. of molecules (stoichiometry)');
    ylabel('Probability density')
    
catch ME1 % if the fit failed (shape can be not gamma distrib-like)
    a_fit = params_gamma(1);
    b_fit = params_gamma(2);
    c_fit = 0;
    rsq_fit = 0;
    stDev_a = (pci(2,1)-pci(1,1))/2;
    stDev_b = (pci(2,2)-pci(1,2))/2;
    stDev_c = 0;
end
    
hold off; % for figure
    
% Save current figure (stoichiometry prob density from all tracks) as .png in 'results_allTracks' folder:
cd(output_folder_path);
figName4 = strcat(short_data_set_label,colour_chosen,'_stoichiom_ProbDensity'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName4)

[mean_gamma,variance_gamma] = gamstat(a_fit,b_fit); % returns mean of and variance for the gamma distribution with parameters specified by params_gamma(1) and params_gamma(2).
% disp('Gamma stats from new fit: ')
% disp(['params gamma distrib:  ',num2str([a_fit b_fit])]);
% disp(['mean_gamma: ',num2str(mean_gamma)]);
% disp(['std_gamma: ',num2str(sqrt(variance_gamma))]);
% disp(['mode gamma: ',num2str((a_fit-1)*b_fit)]);

% numbers for output:
allTrack_stoichiom.a_fit = a_fit;
allTrack_stoichiom.error_a_fit = stDev_a;
allTrack_stoichiom.b_fit = b_fit;
allTrack_stoichiom.error_b_fit = stDev_b;
allTrack_stoichiom.c_fit = c_fit;
allTrack_stoichiom.error_c_fit = stDev_c;
allTrack_stoichiom.rsq_fit = rsq_fit;
allTrack_stoichiom.mean_gamma = round(mean_gamma);
allTrack_stoichiom.mode_gamma_fit = round((a_fit-1)*b_fit); % stoichiom value (x-value) corresponding to the maximum of the fitted curve.
allTrack_stoichiom.std_gamma = round(sqrt(variance_gamma));
allTrack_stoichiom.mode_data = round(vector_stoichiom(stoich_prob_density==max(stoich_prob_density))); % Stoichiom value (x-value) corresponding to the maximum of the data curve.


% ------------------------
% ------------------------
% SAVE INPUTS and info to output variable:
allTrack_inputs.data_set_label = short_data_set_label;
allTrack_inputs.folder_path = initial_folder_path;
% ---------------------
% Use when image labels are only numbers:
% allTrack_inputs.image_numbers = mat2str(image_numbers_1); % convert to a string.
% ---------------------
% % Use when there are image numbers such as '02':
allTrack_inputs.image_numbers = image_numbers_1; % convert to a string.
% ---------------------
allTrack_inputs.minNumPointsInTrack = minNumPointsInTrack;
allTrack_inputs.max_framesAwayFromTimeOrigin = max_framesAwayFromTimeOrigin;
allTrack_inputs.colour_chosen = colour_chosen;
allTrack_inputs.tau_fixed = tau_fixed;
allTrack_inputs.min_rsq_fit_I = min_rsq_fit_I;
allTrack_inputs.max_fit_error = max_fit_error;
allTrack_inputs.prob_accept_fit = prob_accept_fit;
allTrack_inputs.I_singleMolec = I_singleMolec;
allTrack_inputs.error_IsingleMolec = error_IsingleMolec;
allTrack_inputs.num_tracks_used = num_tracks_used; % final number of accepted tracks.

% ----------------------------------
% ----------------------------------
% NOISE ANALYSIS:
% CHECKING NOISE LEVELS FOR DATA SET:

N1 = allTrack_results1.NumDataPoints;   % column vector, no. of data points in track.
sigmaI1 = allTrack_results1.all_I_approx_errorBar; % column vector, standard deviation of flattened I trace (subtracting fit to exp with offset from data).
Iavg1 = allTrack_results1.all_Iavg; % column vector, mean intensity of all points in track (used later on for short tracks).

% Use only tracks with less than 26 points (short enough tracks so that the
% Iavg value is meaningful), so less than 1 second long for oxphos data:
max_Npoints_noiseAnalysis = 25; % default: (25)
Iavg2 = Iavg1(N1<=max_Npoints_noiseAnalysis);
sigmaI2 = sigmaI1(N1<=max_Npoints_noiseAnalysis);
sigmaI2_sq = sigmaI2.^2; % sigma noise squared

% Fit sigma noise squared, sigmaI2_sq versus Iavg2 to a line:
fun_to_fit = fittype('a + b*x','independent','x'); % define exponential funtion to fit to, with 'x' as independent variable.
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
options.StartPoint = [0 max(sigmaI2_sq)/max(Iavg2)]; % give guess parameters for fit. This avoids a warning message. Give in right order!. 
% options.Lower = [0 1]; % give lower bound for parameters. 
% options.Upper = [0 1]; % give upper bound for parameters. 
[fit_result gof] = fit(Iavg2,sigmaI2_sq,fun_to_fit,options); % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% fit_param_names = coeffnames(fit_result); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
intercept_noiseFit = fit_param_values(1); % intercept, bgnd sigma noise squared.
slope_noiseFit = fit_param_values(2); % slope from fit.
rsq_noiseFit = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_intercept_noiseFit = errorSTDEV(1);
stDev_slope_noiseFit = errorSTDEV(2);

% Calculate scaling factor alpha:
alpha_noiseFit = error_IsingleMolec^2/(I_singleMolec*slope_noiseFit);
% New, re-scaled single-molec intensity value:
I_singleMolec_new = I_singleMolec/alpha_noiseFit;


% Save numbers in variable for output:
allTrack_noiseAnalysis.intercept_noiseFit = intercept_noiseFit;
allTrack_noiseAnalysis.stDev_intercept_noiseFit = stDev_intercept_noiseFit;
allTrack_noiseAnalysis.slope_noiseFit = slope_noiseFit;
allTrack_noiseAnalysis.stDev_slope_noiseFit = stDev_slope_noiseFit;
allTrack_noiseAnalysis.rsq_noiseFit = rsq_noiseFit;
allTrack_noiseAnalysis.alpha_noiseFit = alpha_noiseFit; % scaling factor.
allTrack_noiseAnalysis.I_singleMolec_new = I_singleMolec_new; % new fluorophore single molecule value, re-scaled by dividing by alpha_noiseFit.
allTrack_noiseAnalysis.new_stoichiom_mode_value = alpha_noiseFit*allTrack_stoichiom.mode_data; % new stoichiometry value, re-scaled by factor alpha_noiseFit.


% Plot again and hold:
plot(Iavg2,sigmaI2_sq,'*k') % curve in blue
hold on; 
plot(fit_result,'-r'); % plot fit as a red curve
xlabel('Iavg for track');
ylabel('noise-level^2 (stdev^2) for track');
title('Noise data analysis, scaling factor')
hold off;

% Save current figure (noise analysis and scaling factor) as .png in 'results_allTracks' folder, and close figure:
cd(output_folder_path);
figName4 = strcat(short_data_set_label,colour_chosen,'_noiseAnalysisScaleFactor'); % name of figure file to save to.
saveFigurePNG(output_folder_name2,figName4)

% --------------------------------------------


%% Output results:
% ===============================================
% Save output variables onto a .mat file in folder 'results_allTracks':
results_name = strcat('results_',short_data_set_label,colour_chosen); % name of output .mat file.
cd(output_folder_name2)
save(results_name,'results','allTrack_inputs','allTrack_results1','allTrack_results2','allTrack_results3','allTrack_stoichiom','allTrack_noiseAnalysis') % Save variables results, allTrack_results1 and allTrack_results2.


%% Save output results to an ouput excel file inside folder 'results_allTracks':

% Save excel file inside folder "output_folder_name2", named 'results_allTracks':
output_filename = strcat(short_data_set_label,colour_chosen,'allTracks_I0','.xls'); % name of excel file to save to.
% % -----------------
% % When the name of the excel file exceeds 218 characters: error. Change
% % 'bottom' to 'bot':
% output_filename = strcat(short_data_set_label,'bot','allTracks_I0','.xls'); % name of excel file to save to.
% % -----------------
warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.

dataForSheet1 = [fieldnames(allTrack_inputs) struct2cell(allTrack_inputs)]; % numbers, first column is fieldnames, second one is values.
dataForSheet2 = [fieldnames(allTrack_results1)'; num2cell(cell2mat(struct2cell(allTrack_results1)'))]; % column vectors
dataForSheet3 = [fieldnames(allTrack_results2)'; num2cell(cell2mat(struct2cell(allTrack_results2)'))]; % column vectors
dataForSheet4 = [fieldnames(allTrack_results3)'; num2cell(cell2mat(struct2cell(allTrack_results3)'))]; % column vectors
dataForSheet5 = [fieldnames(allTrack_stoichiom) struct2cell(allTrack_stoichiom)]; % numbers, first column is fieldnames, second one is values.
dataForSheet6 = [fieldnames(allTrack_I0) struct2cell(allTrack_I0)]; % numbers, first column is fieldnames, second one is values.
dataForSheet7 = [fieldnames(allTrack_noiseAnalysis) struct2cell(allTrack_noiseAnalysis)]; % numbers, first column is fieldnames, second one is values.

% The order here sets the order of the sheets in the excel file:
xlswrite(output_filename,dataForSheet1,'Inputs'); % write data to sheet 'Inputs' in excel file.
xlswrite(output_filename,dataForSheet2,'Data all Tracks'); % write data to sheet 'Data all Tracks' in excel file.
xlswrite(output_filename,dataForSheet3,'I0 prob density'); % write data to sheet 'I0 prob density' in excel file.
xlswrite(output_filename,dataForSheet6,'Avg I0'); % write data to sheet 'Avg I0' in excel file. Final I0 value obtained from all tracks in data set.
xlswrite(output_filename,dataForSheet4,'stoichiom prob density'); % write data to sheet 'stoichiom prob density' in excel file.
xlswrite(output_filename,dataForSheet5,'Avg stoichiom'); % write data to sheet 'Avg stoichiom' in excel file. Final stoichiometry value obtained from all tracks in data set.
xlswrite(output_filename,dataForSheet7,'Noise analysis, scaling factor'); % write data to sheet 'Noise analysis, scaling factor' in excel file. Final stoichiometry value obtained from all tracks in data set.

% -----------------------------
% Go back to initial directory:
cd(initial_folder_path); 
% cd('..'); % go back to previous directory.
