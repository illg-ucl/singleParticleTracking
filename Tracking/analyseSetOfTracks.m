function  analyseSetOfTracks(minNumPointsInTrack,minNumToAverage,fit_with_weights,colour_chosen)
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Analyse a set of tracks by opening the corresponding excel files and
% building an average track (intensity vs time) with the info in them.
% Then fit the avg track to exponential curves and plot results and display
% fit results.
%
% Inputs:
% - minNumPointsInTrack: min. no. of points a track must have in order to be
% used to form an average track. Default value is 25.
% - minNumToAverage: min no. of tracks which are averaged in order to obtain
% a point in the average track. Default value is 5-10. This value should be
% at least 2 or program will fail.
% - fit_with_weights: parameter which determines if weights are used for
% fits (1) or not (0).
% - colour_chosen: this is a string. It has to be either 'top' or bottom'.
% Choose if you want to look at Green or Red spots (for dual label strains):
% The bottom half of the image corresponds to the Green channel -> 'bottom'.
% The top half of the image corresponds to the Red channel -> 'top'.
% 
% Example of how to run this function:
% analyseSetOfTracks(25,2,1,'top');
% analyseSetOfTracks(25,5,0,'bottom');

% Error control:
if minNumToAverage <2
    disp('The second input needs to be >=2... Change input and re-run function.');
    return
end


%% Build names of folders in which to find tracks to analyse:

% SINGLE LABEL STRAINS:
% ---------------------
% data_set_label_1 = 'ATPase-GFP_'; % string that labels the data set and is part of the name of folders.
% % image_numbers_1 = [89];
% image_numbers_1 = [89 97 99 110 112 116 122 126 128 130 134 136];
 
% data_set_label_1 = 'GFP-nuoF_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [154 156 158 160 162 164 166 168 170  174];

% data_set_label_1 = 'cydB-mCherry_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [228 249 253 257 261 263 265];

% data_set_label_1 = 'mCherry-sdhC_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [279 281 283 285 289 291 293 295 301 303];

% data_set_label_1 = 'cyoA-mCherry_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [309 313 315 317 320 322 324 326];

% data_set_label_1 = 'mCherry-nuoF_'; % string that labels the data set and is part of the name of folders
% image_numbers_1 = [347 349 351 353 357 359];

% DUAL LABEL STRAINS:
% -------------------
% data_set_label_1 = 'cybD-mCherry-ATPase-GFp_'; % string that labels the data set and is part of the name of folders.
% image_numbers_1 = [498 500 509 513 515 518 522 524];

% data_set_label_1 = 'GFP-nuoF-mCherry-sdhC_'; 
% image_numbers_1 = [533 545 547 549 551 553];

% data_set_label_1 = 'cydB-mCherry-GFPuv4-nuoF_'; 
% image_numbers_1 = [427 429 431 438 442 444];

data_set_label_1 = 'cydB-mCherry-GFPuv4-nuoF_'; 
image_numbers_1 = [472 474 478 484 486];


folder_names_1 = cell(1,length(image_numbers_1)); % initialise empty cell array (row).
for i=1:length(image_numbers_1)
    folder_name = strcat(data_set_label_1,num2str(image_numbers_1(i)));
    folder_names_1{1,i} = folder_name;
end

% data_set_label_2 = 'cydB-GFP_'; % string that labels the data set and is part of the name of folders.
% % image_numbers_2 = [02];
% image_numbers_2 = {'02' '11' '13' '15' '25' '29' '37'};
% 
% folder_names_2 = cell(1,length(image_numbers_2)); % initialise empty cell array (row).
% for i=1:length(image_numbers_2)
%     folder_name = strcat(data_set_label_2,image_numbers_2{i});
%     folder_names_2{1,i} = folder_name;
% end

all_folder_names = folder_names_1;
% all_folder_names = folder_names_2; 


%% Go through track files (xls files) in each of the previous folders:

% Initialise index of accepted tracks to accummulate for further analysis:
m = 1; % This is the number of tracks with at least minNumPointsInTrack points in them.

for j =1:length(all_folder_names)
    
    cd(all_folder_names{j}); % move into folder (directory);
    
    xlsFileNames0 = dir('*.xls'); % Get names of excel files in that folder (track analysis .xls files).
    % Error control:
    if isempty(xlsFileNames0) % If there are no .xls track data files, show error and exit function:
        error('Check you are in the correct directory. No .xls files found in folder.');
    end
    
    xlsFileNames = {xlsFileNames0.name}; % cell array of strings.
    
    % Loop through each track analysis xls file:
    for k=1:length(xlsFileNames)
        
        disp(xlsFileNames{k});
        
        % Open and read the data:
        
        % Import the data in the sheet named 'Track info':
        [numeric,txt,raw] = xlsread(xlsFileNames{k},'Track info');
        % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
        str_TrackInfo = cell2struct(raw(:,2),raw(:,1),1);
        NumDataPoints = str_TrackInfo.NumDataPoints;
        AbsTimeOrigin = str_TrackInfo.AbsTimeOrigin; % time origin for a given data set (frame for which laser is switched on).
        TrajStartTime = str_TrackInfo.TrajStartTime; % time at the start of a particular track.
        TimeBetweenFrames = str_TrackInfo.TimeBetweenFrames; % in seconds.
        TopOrBottom = str_TrackInfo.TopOrBottom; % 'top' or 'bottom' region of image (colour channel, red or green);
        % track_with_jumps_flag = str_TrackInfo.track_with_jumps_flag;
        
        
        % Use a criterium to accept which track data we accummulate to use it later:
        if NumDataPoints >= minNumPointsInTrack && ... % if num of points in track is at least minNumPointsInTrack (25):
                strcmp(TopOrBottom,colour_chosen) == 1 % if the track is of the right colour/channel
            
            % Import the data in the sheet named 'Results I fits':
            [numeric,txt,raw] = xlsread(xlsFileNames{k},'Results I fits');
            % Turn imported data from excel file into a structure:
            str_ResultsIfits = cell2struct(raw(:,2),raw(:,1),1);
            I0_fit = str_ResultsIfits.I0_fit; % from fit with no offset
            tau_fit = str_ResultsIfits.tau_fit;
            % Calculate intensity corresponding to first point in track (TrajStartTime),
            % from the fit to an exponential with no offset:
            I_atTrackStart = I0_fit*exp(-(TrajStartTime-AbsTimeOrigin)/tau_fit);
            
            
            % Import the data in the sheet named 'cell coordinates':
            [numeric,txt,raw] = xlsread(xlsFileNames{k},'cell coordinates');
            % Turn imported data from excel file into a structure:
            str_cellCoords = cell2struct(raw(:,2),raw(:,1),1);
            traj_in_pole = str_cellCoords.traj_in_pole;
            cell_width_nm = str_cellCoords.cell_width_nm;
            cell_length_nm = str_cellCoords.cell_length_nm;
            
            % Import the data in the sheet named 'Track data':
            [numeric,txt,raw] = xlsread(xlsFileNames{k},'Track data');
            % Turn imported data from excel file into a structure:
            str_TrackData = struct; % create empty structure to fill up.
            for i = 1:length(txt)
                str_TrackData = setfield(str_TrackData, txt{i}, numeric(:,i)); % create field in the structure.
                % Each field in the structure contains a column vector with the data.
            end
            
            
            % Vectors useful for later:
            % Time in seconds relative to the start of the track:
            time_rel_to_track_start = str_TrackData.timeabs-TrajStartTime; % column vector
            frame_rel_to_track_start = time_rel_to_track_start./TimeBetweenFrames+1; % column vector with frame numbers (there can be jumps of one frame within!).
            % integrated intensity normalised to the intensity at the start of the track:
            Inormalised_to_track_start = str_TrackData.intensity./I_atTrackStart; % column vector
            
            % Accummulate selected data in a matrix:
            for i = 1:length(frame_rel_to_track_start)
                row = int16(frame_rel_to_track_start(i));
                warning off 'MATLAB:intConvertNonIntVal'  % Turn this warning off: Warning: Conversion rounded non-integer floating point value to nearest int16 value.
                column = m;
                % fill up a matrix:
                Inormalised_all(row,column) = Inormalised_to_track_start(i);
                % Matrix "Inormalised_all" contains several columns, each
                % of which is a normalised intensity trace (for one excel file), with the time starting at zero.
                % The row number corresponds to the frame number.
                % Positions with no data (time jumps) appear as zeros (holes).
            end
            
            m = m+1; % advance index for accummulating track info.
            
        end % if NumDataPoints >= minNumPointsInTrack.
        
    end
    
    cd('..'); % go back to previous directory.
    
end

disp(' '); % empty line.
disp(['The number of tracks with at least ',num2str(minNumPointsInTrack),' points is: ',num2str(m-1),'']);


%% Produce a normalised and average intensity track and use it to find an average photobleaching time of a given fluorescent protein:

n =1; % Initialise index for position in final result
for i = 1:size(Inormalised_all,1) % loop through rows in matrix of intensity traces Inormalised_all:
    % the row number, i, corresponds to frame number.
    Inormalised_row0 = Inormalised_all(i,:); % take that row of different intensity values for a given relative time.
    Inormalised_row = Inormalised_row0(Inormalised_row0 ~= 0); % exclude values equal to zero (empty data positions).
    if length(Inormalised_row) >= minNumToAverage % if there are at least 4 elements different from zero in that row:
        Inormalised_avg(n) = mean(Inormalised_row); % mean intensity value for that row, ie, for that time (relative time).
        Inormalised_error(n) = std(Inormalised_row)/sqrt(length(Inormalised_row)); % use standard deviation/sqrt(num points) as standard error of the mean.
        % Inormalised_error(n) = std(Inormalised_row); % use standard deviation as error of the mean.
        time_axis(n) = (i-1)*TimeBetweenFrames; % the row number in matrix Inormalised_row, which can have holes, is the frame number (i).
        n = n+1;
        % Note that there can be jumps in time on the averaged-normalised intensity trace!
    end
    
end

% Error control:
if exist('Inormalised_avg','var')==0 % if the variable does not exist
    disp('Not enough tracks to average over. Try to change input parameters.');
    return % exit program
end

% Data for fits. We fit the averaged and normalised integrated spot intensity (bgnd subtracted) to an exponential:
IforFit = Inormalised_avg'; % column vector
tforFit = time_axis'; %  column vector,  time relative to start of track.
% Weights for fit:
if fit_with_weights == 1
    weights_forFit = (1./Inormalised_error).^2./sum((1./Inormalised_error).^2); % normalised weights for fit (row vector (fit works like this)).
else
    weights_forFit = []; % for no weighing in fits, make this an empty vector.
end


%% Fit intensity data to an exponentially decaying function (with no offset):

exp_no_offset = fittype('I0*exp(-t/tau)','independent','t'); % define exponential funtion to fit to, with 't' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = 1; % since I has been normalised.
guess_tau = 2; % in seconds;
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 0]; % Lower bounds for fit parameters. In order: I0, tau.
options.Upper = [2 30]; % Upper bounds for fit parameters. In order: I0, tau.
options.Weights = weights_forFit;
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
options.Weights = weights_forFit;
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

errorbar(time_axis,Inormalised_avg,Inormalised_error,'.k'); % plot with error bars.
hold on;
plot(fit_result_I,'b'); % plot exponential fit (no offset) as BLUE line.
plot(fit_result_I_wo,'r'); % plot exponential fit (with offset) as RED line.
legend('avgd data','exp fit-no offset','exp fit-with offset');
% legend('hide');
ylim([0 1.2]); 
xlim([0 n*TimeBetweenFrames]);
xlabel('time from track start (s)'); 
ylabel('Averaged normalised intensity');
hold off;



