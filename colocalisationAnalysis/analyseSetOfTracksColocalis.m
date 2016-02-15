function  analyseSetOfTracksColocalis(minNumPointsInTrack,minNumFramesOverlap)
%
% Created by Isabel Llorente-Garcia, July 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Look for colocalised tracks. Find tracks with data points at the same
% times and in different channels (green/bottom and red/top).
%
% Analyse a set of tracks by opening the corresponding excel files in a
% chosen folder (for a given data set and a given image sequence) and using the info in them.
% Use all tracks with a number of points at least equal to
% minNumPointsInTrack. Go track by track.
% Do this for all the folders (one per image sequence) corresponding to a
% given data set.
% This function only makes sense for dual label strains.
%
% -- Inputs:
% - minNumPointsInTrack: min. no. of points a track must have in order to be
% considered. (Default value is 5.)
% - minNumFramesOverlap: min no. of frames (data points) overlapping for
% two tracks of different channels to be analysed for co-localisation (default is 2).
%
% -- Outputs: 
% One excel file and png plot is produced per time-overlapping track pair.
% The excel file contains Sheets with the values from structure "result_numbers".
% First element is channel j=1 (top/red), second element is channel j = 2 (bottom/green).
% Other sheets contain the "result_vectors" structure, also for j=1 and
% j=2. 
%
% Example of how to run this function:
% analyseSetOfTracksColocalis(5,2);


%% Build names of folders in which to find tracks to analyse:

% For building paths in which to find track excel files:
% NOTE to user: uncomment one of the following options: 

% DUAL LABEL STRAINS:
% -------------------
% short_data_set_label = 'cybD-mCherry-ATPase-GFp_'; % string that labels the data set and is part of the name of folders.
% top_folder = 'cybD-mCherry-ATPase-GFp_I0_top_440pm250'; % folder name (within data set) which contains top channel accepted trajectories.
% bottom_folder = 'cybD-mCherry-ATPase-GFp_I0_bottom_1000pm700'; % folder name (within data set) which contains bottom channel accepted trajectories.
% image_labels = {'498' '500' '509' '513' '515' '518' '522' '524'};
% % image_labels = {'498'}; 
% Displacement between top and bottom channels, in order of images, [ydisplacement xdisplacement] = overlapBF(...), from the corresponding bright field image, or by hand if that is no good.
% Top channel shifted by:
% translateBy_toUse = {[5 -4] [6 -5] [4 -7] [6 -1] [3 -5] [2 -5] [3 -6] [5 -5]}; 
% % translateBy_toUse = {[5 -4]}; % use this only for image_labels = {'498'}; 
% ---------------
% short_data_set_label = 'GFP-nuoF-mCherry-sdhC_'; % string that labels the data set and is part of the name of folders.
% top_folder = 'GFP-nuoF-mCherry-sdhC_I0_top_440pm250'; % folder name (within data set) which contains top channel accepted trajectories.
% bottom_folder = 'GFP-nuoF-mCherry-sdhC_I0_bottom_1000pm700'; % folder name (within data set) which contains bottom channel accepted trajectories.
% image_labels = {'533' '545' '547' '549' '551' '553'};
% Displacement between top and bottom channels, in order of images, [ydisplacement xdisplacement] = overlapBF(...), from the corresponding bright field image, or by hand if that is no good.
% Top channel shifted by:
% translateBy_toUse = {[2 -5] [0 -5] [-2 -4] [ -1 -3] [3 -4] [0 1]}; 
% ---------------
% short_data_set_label = 'cydB-mCherry-GFPuv4-nuoF_'; 
% top_folder = 'cydB-mCherry-GFPuv4-nuoF_I0_top'; % folder name (within data set) which contains top channel accepted trajectories.
% bottom_folder = 'cydB-mCherry-GFPuv4-nuoF_I0_bot'; % folder name (within data set) which contains bottom channel accepted trajectories.
% image_labels = {'427' '429' '431' '438' '442' '444'};
% Displacement between top and bottom channels, in order of images, [ydisplacement xdisplacement] = overlapBF(...), from the corresponding bright field image, or by hand if that is no good.
% Top channel shifted by:
% translateBy_toUse = {}; 
% ---------------
short_data_set_label = 'cydB-mCherry-GFPuv4-nuoF_'; 
top_folder = 'cydB-mCherry-GFPuv4-nuoF_I0_top'; % folder name (within data set) which contains top channel accepted trajectories.
bottom_folder = 'cydB-mCherry-GFPuv4-nuoF_I0_bot'; % folder name (within data set) which contains bottom channel accepted trajectories.
image_labels = {'472' '474' '478' '484' '486'};
% Displacement between top and bottom channels, in order of images, [ydisplacement xdisplacement] = overlapBF(...), from the corresponding bright field image, or by hand if that is no good.
% Top channel shifted by:
translateBy_toUse = {[6 -4] [4 -4] [6 -4] [6 -7] [5 -6]}; 
% ---------------

all_folder_names = {top_folder,bottom_folder};


%% Create new directory for saving plots and results:

% Make new folder (new directory) to save trajectory analysis results:
initial_folder_path = cd;
% Create new folder for outputs inside current folder:
output_folder_name = ('colocResultsTrackByTrack'); % Colocalisation results.
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

disp(' ') % empty line.
disp('------------------------------')
disp(['Data set: ' short_data_set_label])

% Loop through all image sequences in a data set:
for i = 1:length(image_labels) % i is the index for each image sequence.
    
    % Translate top channel by:
    shiftTopChannelBy_y = translateBy_toUse{i}(1);
    shiftTopChannelBy_x = translateBy_toUse{i}(2);
    shiftBottomChannelBy_y = 0; % Do not translate bottom channel.
    shiftBottomChannelBy_x = 0; % Do not translate bottom channel.
    
    disp('--------')
    disp(['Image sequence: ' image_labels{i}])
    
    % Go through folders for both top and bottom accepted tracks:
    for j = 1:length(all_folder_names) % j=1,2 is the index for each channel in an image, ie., top/red and bottom/green.
        
        cd(all_folder_names{j}); % move into folder (directory) containing all top or bottom accepted tracks.
       
        xlsFileNames0 = dir(strcat('*',image_labels{i},'*.xls')); % Get names of excel files in that folder (accepted track-analysis .xls files).
        
        % ------------
        % Initialise index of accepted tracks to accummulate for further analysis:
        m = 1;
        % This is the number of tracks with at least
        % minNumPointsInTrack points in them.
        
        % Error control:
        if ~isempty(xlsFileNames0) % If there are some .xls track-data files for that image sequence:

            xlsFileNames = {xlsFileNames0.name}; % cell array of strings.
            
            % Loop through each track-analysis xls-file inside folder j:
            for k=1:length(xlsFileNames) % k is the index for each excel file.
                
                % disp(xlsFileNames{k});
                
                % pos1 = strfind(xlsFileNames{k},'_bis.xls'); % position of the start of the string '_bis.xls' in the excel file name.
                % figName_0 = xlsFileNames{k}(1:(pos1-1)); % Take the name of the input excel file (with the end bit '_bis.xls' removed) as the new name for saving a later .png figure.
                % figName = strcat(figName_0,'_coloc'); % name of figure file to save to.
                
                % Open and read the data:
                % -----------------------
                % Import the data in the sheet named 'Track results':
                [numeric,txt,raw] = xlsread(xlsFileNames{k},'Track info');
                % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
                str_TrackInfo = cell2struct(raw(:,2),raw(:,1),1);
                NumDataPoints = str_TrackInfo.NumDataPoints;
                
                if NumDataPoints >= minNumPointsInTrack % Use track only if it has at least "minNumPointsInTrack" data points.
                    
                    % Import the data in the sheet named 'Track data':
                    [numeric1,txt1,raw1] = xlsread(xlsFileNames{k},'Track data');
                    % Turn imported data from excel file into a structure:
                    str_TrackData = struct; % create empty structure to fill up.
                    for u = 1:length(txt1)
                        str_TrackData = setfield(str_TrackData, txt1{u}, numeric1(:,u)); % create field in the structure.
                        % Each field in the structure contains a column vector with the data.
                    end
                    
                    % Save relevant track info in a structure to later find
                    % tracks which occur at same time in different channels (j):
                    s(j,m).xlsFileName = xlsFileNames{k};
                    s(j,m).FirstTrajFrame = str_TrackInfo.FirstTrajFrame; % First frame in this particular track.
                    s(j,m).LastTrajFrame = str_TrackInfo.LastTrajFrame; % First frame in this particular track.
                    s(j,m).TopOrBottom = str_TrackInfo.TopOrBottom; % 'top' or 'bottom' region of image (colour channel, red or green);
                    s(j,m).NumDataPoints = NumDataPoints;
                    s(j,m).frames_vector = str_TrackData.frame; % column vector with frame numbers for all points in track.
                    
                    m = m+1;
                end
                
            end % end of k loop (each track).
            
            m_final(i,j) = m - 1; % final total number of long enough tracks for each image i and channel j.
                        
        end % end of "if" there are tracks for that image sequence.
        
        
        % cd('..'); % go back to previous directory.
        cd(initial_folder_path) % go to initial directory where folder has been created to save output.
        
    end % end of j loop (top/red or bottom/green channels).
    
    disp(['There are ' num2str(m_final(i,1)) ' tracks with at least ' num2str(minNumPointsInTrack) ' points from top/red channel.'])
    disp(['There are ' num2str(m_final(i,2)) ' tracks with at least ' num2str(minNumPointsInTrack) ' points from bottom/green channel.'])
    
    % ------------
    % Error control: if for that video there are only tracks of one colour:
    if m_final(i,1)==0 || m_final(i,2)==0
        disp(['For image sequence ',image_labels{i},', tracks of only one colour (one channel) were found.'])
        disp('Colocalisation analysis is meaningless. Skipped.')
    else % if there are tracks of more than one colour (main case):
        
        % Loop through long enough tracks and find the ones from different
        % channels that overlap in time, save results to variable "r":
        n = 1; % initialise index of track-pairs with some overlap in time:
        for k1 = 1:m_final(i,1) % loop through long enough tracks from first channel (top).
            for k2 = 1:m_final(i,2) % loop through long enough tracks from second channel (bottom).
                % For top channel:
                vector_frames_top = s(1,k1).frames_vector; % vector of frame numbers in track, with no jumps, top channel.
                % For bottom channel:
                vector_frames_bot = s(2,k2).frames_vector; % vector of frame numbers in track, with no jumps, bottom channel.
                % Find intersection of the frame-number vectors for the different channels:
                [c top_pos bot_pos] = intersect(vector_frames_top,vector_frames_bot);
                if length(c) >= minNumFramesOverlap % if intersection has at least "minNumFramesOverlap" elements, ie., if there is long enough time overlap between the different channel tracks:
                    % Save file names of time overlapping, top and bottom channel tracks
                    r(n).filename_top = s(1,k1).xlsFileName;
                    r(n).filename_bot = s(2,k2).xlsFileName;
                    r(n).frames_top = vector_frames_top; % column vector
                    r(n).frames_bot = vector_frames_bot; % column vector
                    r(n).overlap_frames = c; % common frame numbers for top and bottom tracks, column vector.
                    r(n).top_pos = top_pos; % position of overlapping frames in vector_frames_top, column vector.
                    r(n).bot_pos = bot_pos; % position of overlapping frames in vector_frames_bot, column vector.
                    n = n + 1; 
                end
            end
        end
        n_final(i) = n - 1; % final number of track-pairs which overlap in time for a given image sequence (i).

        disp(['There are ' num2str(n_final(i)) ' track-pairs which overlap in time for image sequence ' image_labels{i} '.'])
        
        % -------------
        % Loop through time-overlapping pairs of tracks:
        if n_final(i) > 0 % if there are are time-overlapping tracks;
            
            for p = 1:n_final(i)   % Loop through time-overlapping track pairs:
                
                % Filenames of the two time-overlapping tracks in the pair:
                overlapping_track_filenames = {r(p).filename_top,r(p).filename_bot};
                position_overlapping_frames = [(r(p).top_pos)'; (r(p).bot_pos)']; % first row is for channel 1 (top), second raw is for channel 2 (bottom).
                % SAVE TO OUTPUT: create a "result_numbers" structure which will be saved to an
                % excel file. First element (result_numbers(1)) is channel 1 (top/red),
                % second element (result_numbers(2)) is channel 2 (bottom/green).
                % There will also be a "result_vectors" structure saved
                % later.
                % One excel file saved per time-overlapping track pair.
                result_numbers(1).filename = r(p).filename_top; 
                result_numbers(1).originalNumFrames = length(r(p).frames_top);
                result_numbers(1).numOverlappingFrames = length(r(p).overlap_frames); % same for both channels.
                result_numbers(2).filename = r(p).filename_bot;
                result_numbers(2).originalNumFrames = length(r(p).frames_bot);
                result_numbers(2).numOverlappingFrames = length(r(p).overlap_frames);
                % --------------------------
                % READ IN MORE DETAILED DATA from the .xls filenames for
                % each track in the pair:
                for j = 1:length(all_folder_names) % j=1,2 is the index for each channel in an image, ie., top/red and bottom/green.
                    cd(all_folder_names{j}); % move into folder (directory) containing all top or bottom accepted tracks.
                    % --------------
                    % Import the data (numbers) in the sheet named 'Track info':
                    [numeric,txt,raw] = xlsread(overlapping_track_filenames{j},'Track info');
                    % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
                    str_TrackInfo = cell2struct(raw(:,2),raw(:,1),1);
                    result_numbers(j).ImageSizeVerticPix = str_TrackInfo.ImageSizeVerticPix; % frame vertical size in pixels
                    result_numbers(j).ImageSizeHorizPix = str_TrackInfo.ImageSizeHorizPix; % frame horizontal size in pixels.
                    result_numbers(j).TopOrBottom = str_TrackInfo.TopOrBottom; % 'top' or 'bottom' region of image (colour channel, red or green);
                    result_numbers(j).FirstTrajFrame = str_TrackInfo.FirstTrajFrame; % First frame in this particular track.
                    result_numbers(j).LastTrajFrame = str_TrackInfo.LastTrajFrame; % First frame in this particular track.
                    result_numbers(j).TimeBetweenFrames = str_TrackInfo.TimeBetweenFrames; % in seconds.
                    result_numbers(j).pixelsize_nm = str_TrackInfo.pixelsize_nm;
                    result_numbers(j).meanSpotWidth = str_TrackInfo.meanSpotWidth;
                    % Save inputs:
                    result_numbers(j).minNumPointsInTrack = minNumPointsInTrack;
                    result_numbers(j).minNumFramesOverlap = minNumFramesOverlap;   
                    if j==1 % amounts by which to shift top channel along x and y:
                        result_numbers(j).shiftChannelBy_y = shiftTopChannelBy_y;
                        result_numbers(j).shiftChannelBy_x = shiftTopChannelBy_x;
                    elseif j==2 % amounts by which to shift bottom channel (0, ie, left unshifted)
                        result_numbers(j).shiftChannelBy_y = shiftBottomChannelBy_y;
                        result_numbers(j).shiftChannelBy_x = shiftBottomChannelBy_x;
                    end
                    % --------------
                    % Import the data (column vectors) in the sheet named 'Track data':
                    [numeric1,txt1,raw1] = xlsread(overlapping_track_filenames{j},'Track data');
                    % Turn imported data from excel file into a structure:
                    str_TrackData = struct; % create empty structure to fill up.
                    for u = 1:length(txt1)
                        str_TrackData = setfield(str_TrackData, txt1{u}, numeric1(:,u)); % create field in the structure.
                        % Each field in the structure contains a column vector with the data.                        
                    end
                    % Import all values as column vectors:
                    frame = str_TrackData.frame; % column vector.
                    timeabs = str_TrackData.timeabs; % column vector.
                    intensity = str_TrackData.intensity; % column vector.
                    xvalues0 = str_TrackData.xvalues0; % column vector.
                    yvalues0 = str_TrackData.yvalues0; % column vector.
                    x_values_cell = str_TrackData.x_values_cell; % column vector.
                    y_values_cell = str_TrackData.y_values_cell; % column vector.
                    sigma_IspotFit = str_TrackData.sigma_IspotFit; % sigma from Gaussian fit to spot.
                    I0_IspotFit = str_TrackData.I0_IspotFit; % amplitude I0 from Gaussian fit to spot.                                                            
                    % Take only the positions in those vectors
                    % corresponding to the frames which overlap in time in
                    % both channels:
                    result_vectors(j).frames = frame(position_overlapping_frames(j,:));
                    result_vectors(j).timeabs = timeabs(position_overlapping_frames(j,:));
                    result_vectors(j).intensity = intensity(position_overlapping_frames(j,:));
                    result_vectors(j).xvalues0 = xvalues0(position_overlapping_frames(j,:));
                    result_vectors(j).yvalues0 = yvalues0(position_overlapping_frames(j,:));
                    % Shifted positions of track (shift top, leave bottom channel the same):
                    if j == 1 % shift top channel values by [shiftTopChannelBy_y shiftTopChannelBy_x]
                        result_vectors(j).xvalues0_shifted = shiftTopChannelBy_x + result_vectors(j).xvalues0;
                        result_vectors(j).yvalues0_shifted = shiftTopChannelBy_y + result_vectors(j).yvalues0;
                    elseif j ==2 % leave bottom channel unshifted except for the y-coordinates, which get shifted by half the vertical frame size to be able to compare with the y-coord from the top channel:
                        result_vectors(j).xvalues0_shifted = shiftTopChannelBy_x + result_vectors(j).xvalues0;
                        result_vectors(j).yvalues0_shifted = shiftTopChannelBy_y + result_vectors(j).yvalues0 - round((result_numbers(j).ImageSizeVerticPix)/2);
                    end
                    result_vectors(j).x_values_cell = x_values_cell(position_overlapping_frames(j,:));
                    result_vectors(j).y_values_cell = y_values_cell(position_overlapping_frames(j,:));
                    result_vectors(j).sigma_IspotFit = sigma_IspotFit(position_overlapping_frames(j,:));
                    result_vectors(j).I0_IspotFit = I0_IspotFit(position_overlapping_frames(j,:));
                    
                    % --------------
                    % Import the data (numbers) in the sheet named 'Mobility results':
                    [numeric2,txt2,raw2] = xlsread(overlapping_track_filenames{j},'Mobility results');
                    % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
                    str_Mobility = cell2struct(raw2(:,2),raw2(:,1),1);
                    result_numbers(j).fit_msd_line_rsq = str_Mobility.fit_msd_line_rsq;
                    result_numbers(j).fit_msd_conf_rsq = str_Mobility.fit_msd_conf_rsq;
                    result_numbers(j).diffusion1D_coeff_micronsSqrdPerSec = str_Mobility.diffusion1D_coeff_micronsSqrdPerSec;
                    result_numbers(j).diffusion1D_coeff_error = str_Mobility.diffusion1D_coeff_error;
                    result_numbers(j).size_1Dconfin_nm = str_Mobility.size_1Dconfin_nm;
                    result_numbers(j).size_1Dconfin_nm_error = str_Mobility.size_1Dconfin_nm_error;
                    result_numbers(j).Brownian_flag = str_Mobility.Brownian_flag;
                    result_numbers(j).Confined_flag = str_Mobility.Confined_flag;
                    result_numbers(j).OtherMobility_flag = str_Mobility.OtherMobility_flag;
                    % --------------
                    % Import the data (numbers) in the sheet named 'cell coordinates':
                    [numeric3,txt3,raw3] = xlsread(overlapping_track_filenames{j},'cell coordinates');
                    % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
                    str_cellCoords = cell2struct(raw3(:,2),raw3(:,1),1);
                    result_numbers(j).cell_width_nm = str_cellCoords.cell_width_nm;
                    result_numbers(j).cell_length_nm = str_cellCoords.cell_length_nm;
                    result_numbers(j).traj_in_pole = str_cellCoords.traj_in_pole;
                    % --------------
                    % Import the data (numbers) in the sheet named 'initial_I data':
                    [numeric4,txt4,raw4] = xlsread(overlapping_track_filenames{j},'initial_I data');
                    % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
                    str_stoichiom = cell2struct(raw4(:,2),raw4(:,1),1);
                    result_numbers(j).stoichiom = str_stoichiom.stoichiom;
                    result_numbers(j).error_stoichiom = str_stoichiom.error_stoichiom;
                    % --------------
                    
                    cd(initial_folder_path) % go to initial directory where folder has been created to save output.
                
                end % end of loop through j=1 and j=2, both channels in the track pair.
                
                % ----------------------------------------
                % Calculate overlap integral (vectors):
                % Use parameters of Gaussian fit to bright spots (use normalised Gaussians of amplitude equal to 1).   
                % Index "i" is the image sequence, index "p" is the number
                % of time-overlapping track pair.
                time_abs{i,p} = result_vectors(1).timeabs; % column vector, same for both channels 1 and 2 (overlapping times).
                sigma_top{i,p} = result_vectors(1).sigma_IspotFit; % column vector.
                sigma_bottom{i,p} = result_vectors(2).sigma_IspotFit; % column vector.
                x_top{i,p} = result_vectors(1).xvalues0_shifted; % column vector, shifted x position in top channel.
                y_top{i,p} = result_vectors(1).yvalues0_shifted; % column vector, shifted y position in top channel.
                x_bottom{i,p} = result_vectors(2).xvalues0_shifted; % column vector, shifted x position in bottom channel.
                y_bottom{i,p} = result_vectors(2).yvalues0_shifted; % column vector, shifted y position in bottom channel.
                delta_x{i,p} = x_top{i,p} - x_bottom{i,p}; % column vector, distance between spot positions along x, for each frame. 
                delta_y{i,p} = y_top{i,p} - y_bottom{i,p}; % column vector, distance between spot positions along y. 
                % Normalisation constant (value of overlap integral for
                % perfect overlap of two Gaussian spots):
                normalise_const{i,p} = pi*sigma_top{i,p}.*sigma_bottom{i,p}; % column vector.
                % Overlap integral: overlap of two 2D Gaussians (integral of their product), column vector:
                overlap{i,p} = pi./(1./(2*sigma_top{i,p}.^2)+1./(2*sigma_bottom{i,p}.^2)).*exp(-(delta_x{i,p}.^2+delta_y{i,p}.^2)./(2*(sigma_top{i,p}.^2+sigma_bottom{i,p}.^2)))./normalise_const{i,p}; % vector.
                % Save results for excel file later (all column vectors):
                % Note that some of the things are saved again even if they are in
                % the excel Sheets for each channel already:
                results_overlap_vectors.time_abs = time_abs{i,p};   
                results_overlap_vectors.sigma_top = sigma_top{i,p};
                results_overlap_vectors.sigma_bottom = sigma_bottom{i,p};
                results_overlap_vectors.x_top = x_top{i,p};
                results_overlap_vectors.y_top = y_top{i,p};
                results_overlap_vectors.x_bottom = x_bottom{i,p};
                results_overlap_vectors.y_bottom = y_bottom{i,p};
                results_overlap_vectors.delta_x = delta_x{i,p};
                results_overlap_vectors.delta_y = delta_y{i,p};
                results_overlap_vectors.normalise_const = normalise_const{i,p};
                results_overlap_vectors.overlap = overlap{i,p};
                % ----------------------------------------
                
                % ----------------------------------------
                % OUTPUT to EXCEL file and .png graph files:
                % ----------------------------------------
                cd(output_folder_path); % move into output folder
                % Save result structures to an excel file for that
                % track pair, inside folder "output_folder_name":
                output_filename = strcat(short_data_set_label,image_labels{i},'_trackPair',num2str(p),'.xls'); % name of excel file to save to.
                warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.
                % Prepare data to save:
                dataForSheet1 = [fieldnames(result_numbers(1)) struct2cell(result_numbers(1))]; % numbers, first column is fieldnames, second one is values.
                dataForSheet2 = [fieldnames(result_numbers(2)) struct2cell(result_numbers(2))]; % numbers, first column is fieldnames, second one is values.
                dataForSheet3 = [fieldnames(result_vectors(1))'; num2cell(cell2mat(struct2cell(result_vectors(1))'))]; % column vectors
                dataForSheet4 = [fieldnames(result_vectors(2))'; num2cell(cell2mat(struct2cell(result_vectors(2))'))]; % column vectors
                dataForSheet5 = [fieldnames(results_overlap_vectors)'; num2cell(cell2mat(struct2cell(results_overlap_vectors)'))]; % column vectors
                % The order here sets the order of the sheets in the excel file:
                xlswrite(output_filename,dataForSheet1,'top channel numbers'); % write data to sheet 'top channel numbers' in excel file.
                xlswrite(output_filename,dataForSheet2,'bottom channel numbers'); % write data to sheet 'bottom channel numbers' in excel file.
                xlswrite(output_filename,dataForSheet3,'top channel vectors'); % write data to sheet 'top channel vectors' in excel file.
                xlswrite(output_filename,dataForSheet4,'bottom channel vectors'); % write data to sheet 'bottom channel vectors' in excel file.
                xlswrite(output_filename,dataForSheet5,'Overlap result vectors'); % write data to sheet 'Overlap result vectors' in excel file.
                % ----------------------------------------
                % Create and save .png file:
                h0 = figure('Tag','Track-pair plots','color','white','units','inches','position',[4 2.7 6 3]); 
                subplot(1,2,1)
                plot(time_abs{i,p},overlap{i,p})
                xlabel('time (s)')
                ylabel('overlap fraction')  
                subplot(1,2,2)
                plot(time_abs{i,p},delta_x{i,p},time_abs{i,p},delta_y{i,p})
                xlabel('time (s)')
                ylabel('delta (pix)')
                legend('delta-x','delta-y','Location','NorthOutside');
                % Save figure as .png in folder created (output_folder_name), and close figure:
                cd(initial_folder_path); % Go back to initial directory.
                figName0 = strcat(short_data_set_label,image_labels{i},'_trackPair',num2str(p),'_overlap'); % name of figure file to save to.
                saveFigurePNG(output_folder_name,figName0)               
                % ----------------------------------------
                % Create and save another .png file:
                h1 = figure('Tag','Track-pair plots','color','white','units','inches','position',[4 2.7 16 9]); 
                subplot(2,3,1)
                plot(result_vectors(1).timeabs,result_vectors(1).xvalues0_shifted,result_vectors(2).timeabs,result_vectors(2).xvalues0_shifted)
                xlabel('time (s)')
                ylabel('x0-shifted (pix)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                subplot(2,3,4)
                plot(result_vectors(1).timeabs,result_vectors(1).yvalues0_shifted,result_vectors(2).timeabs,result_vectors(2).yvalues0_shifted)
                xlabel('time (s)')
                ylabel('y0-shifted (pix)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                subplot(2,3,2)
                plot(result_vectors(1).timeabs,result_vectors(1).x_values_cell,result_vectors(2).timeabs,result_vectors(2).x_values_cell)
                xlabel('time (s)')
                ylabel('x-cell (pix)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                subplot(2,3,5)
                plot(result_vectors(1).timeabs,result_vectors(1).y_values_cell,result_vectors(2).timeabs,result_vectors(2).y_values_cell)
                xlabel('time (s)')
                ylabel('y-cell (pix)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                subplot(2,3,3)
                plot(result_vectors(1).xvalues0_shifted,result_vectors(1).yvalues0_shifted,result_vectors(2).xvalues0_shifted,result_vectors(2).yvalues0_shifted)
                xlabel('x0-shifted (pix)')
                ylabel('y0-shifted (pix)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                subplot(2,3,6)
                plot(result_vectors(1).timeabs,result_vectors(1).intensity,result_vectors(2).timeabs,result_vectors(2).intensity)
                xlabel('time (s)')
                ylabel('intensity (arb)')
                legend('top/red channel','bottom/green channel','Location','NorthOutside');
                % ----------------------------------------
                cd(initial_folder_path); % Go back to initial directory.
                % Save figure as .png in folder created (output_folder_name), and close figure:
                figName = strcat(short_data_set_label,image_labels{i},'_trackPair',num2str(p)); % name of figure file to save to.
                saveFigurePNG(output_folder_name,figName)
                
            end % end of loop through time-overlapping track pairs (p).
            
        end % enf of if there are some overlapping track pairs.
        
    end % if there are tracks of more than one colour.
   
    
end % end of i loop (each image sequence in a data set).


% ---------------------------------------------------------------------
% PUT RESULTS FOR ALL TRACK-PAIRS AND FOR ALL IMAGE SEQUENCES TOGETHER:
% ---------------------------------------------------------------------
% The cell array "overlap" contains in each row the results for a given
% image sequence. Each column corresponds to a given track-pair in that
% image sequence. Note: some elements might be empty, the number of columns
% in all rows is the maximun no. of time-overlapping track pairs from those
% for all image sequences. Each element, eg. overlap{1,1} is a column
% vector with all the normalised overlap integral values for a given track
% pair.
overlap_all_values = []; % Initialise empty vector to accummulate all overlap values.
delta_x_all_values = []; % Initialise empty vector to accummulate all overlap values.
delta_y_all_values = []; % Initialise empty vector to accummulate all overlap values.
sigma_ratio_all_values = []; % Initialise empty vector to accummulate all overlap values.
for ii=1:size(overlap,1) % loop through the number of image sequences analised
   for pp=1:size(overlap,2)
       overlap_all_values = [overlap_all_values ; overlap{ii,pp}]; % column vector
       delta_x_all_values = [delta_x_all_values ; delta_x{ii,pp}]; % column vector
       delta_y_all_values = [delta_y_all_values ; delta_y{ii,pp}]; % column vector
       sigma_ratio_all_values = [sigma_ratio_all_values ; sigma_top{ii,pp}./sigma_bottom{ii,pp}]; % column vector
   end
end
results_all_vectors.overlap_all_values = overlap_all_values; % column vector
results_all_vectors.delta_x_all_values = delta_x_all_values; % column vector
results_all_vectors.delta_y_all_values = delta_y_all_values; % column vector
results_all_vectors.sigma_ratio_all_values = sigma_ratio_all_values; % column vector

% ----------------------------------------
% Plot overlap results for all track-pairs and for all image sequences
% together:
% ----------------------------------------
h_all = figure('Tag','Track-pair plots','color','white','units','inches','position',[4 2.7 8 8]); 
subplot(2,2,1)
hist(overlap_all_values,50)
xlabel('overlap fraction all')
ylabel('frequency')
subplot(2,2,2)
hist(sigma_ratio_all_values,25)
xlabel('sigma ratio all')
ylabel('frequency')
subplot(2,2,3)
hist(delta_x_all_values,25)
xlabel('delta-x all')
ylabel('frequency')
subplot(2,2,4)
hist(delta_y_all_values,25)
xlabel('delta-y all')
ylabel('frequency')
% Save figure as .png in folder created (output_folder_name), and close figure:
cd(initial_folder_path); % Go back to initial directory.
figName_all = strcat(short_data_set_label,'_allImages_allTrackPairs'); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName_all)

% -----------------------------------------
% OUTPUT to EXCEL file and .png graph files:
% -----------------------------------------
cd(output_folder_path); % move into output folder
% Save result structures to an excel file for that
% track pair, inside folder "output_folder_name":
output_filename = strcat(short_data_set_label,'_allImages_allTrackPairs','.xls'); % name of excel file to save to.
warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.
% Prepare data to save:
dataForSheet = [fieldnames(results_all_vectors)'; num2cell(cell2mat(struct2cell(results_all_vectors)'))]; % column vectors
% The order here sets the order of the sheets in the excel file:
xlswrite(output_filename,dataForSheet,'results_allImages_allTrackPairs'); % write data to sheet 'top channel vectors' in excel file.

cd(initial_folder_path); % Go back to initial directory.