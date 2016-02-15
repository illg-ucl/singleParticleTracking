function saveVideoManyTracks(image_label,folderWithTrajs,list_trajNumbers,roi_top,roi_bottom,forced_start_frame,forced_end_frame,saveAvi,output_label)
% 
% Created by Isabel Llorente-Garcia, August 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
%
% Create a single video with the image sequence and several tracks overlaid
% on top as they happen. The video goes from the first frame to the last
% which appear in all tracks together.
% 
% INPUTS:
% - image_label: string which labels the image sequence, eg, '533'.
% Note that the current directory must be that where that image
% sequence is.
% - folderWithTrajs: string, name of folder containing the excel files for
% the trajectories/tracks to plot. There is one excel file per trajectory
% in the input folder folderWithTrajs, with names such as eg. 'GFP-nuoF-mCherry-sdhC_533_traj1.xls'.
% - list_trajNumbers: row vector with the numbers of all the tracks we want
% to plot overlaid on the image sequence.
% - roi_top and roi_bottom: regions of interest to use for the top and
% bottom half of the images. Note the roi's have to be of same size.s
% Use 'full_image' if you want to use the whole
% image.
% - forced_start_frame: force the start of the video a particular frame
% number. If 'automatic', then the first frame is the very first frame
% found in all trajectories to plot.
% - forced_end_frame: force the end of the video to correspond to a
% particular frame. If 'automatic', then the last frame is the latest frame
% found in all trajectories to plot.
% - saveVideoAvi: use 1 for yes and 0 for no. If 0, the video will be
% displayed but not shown as an avi file.
% - output_label: string to label the output .avi video file as you wish.
%
% Examples of how to call this function:
% For image 533 of data set GFP-nuoF & mCherry-sdh
% folder1 = 'C:\Isabel\ExperimData\HeikoData\GFP-nuoF & mCherry-sdhC\GFP-nuoF & mCherry-sdhC TIRF\good\GFP-nuoF-mCherry-sdhC_533'; 
% (list1 = [1 3 4 7];)
% folder1 = 'Z:\ExperimData\HeikoData\GFP-nuoF & mCherry-sdhC\GFP-nuoF & mCherry-sdhC TIRF\good\GFP-nuoF-mCherry-sdhC_533'; 
% or  list1 = [1 3 4 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35];
% saveVideoManyTracks('533',folder1,list1,'full_image','full_image','automatic','automatic',0,'');
% showing full image and not saving, no output extra label.
% ----------
% With top and bottom regions of interest, for cell on the left:
% roi_topA = [100 145 52 95];  roi_bottomA = [93 138 314 357];
% saveVideoManyTracks('533',folder1,list1,roi_topA,roi_bottomA,'automatic','automatic',0,'_A');
% ----------
% For cell on the right:
% roi_topB = [406 448 34 70];  roi_bottomB = [405 447 293 329];
% saveVideoManyTracks('533',folder1,list1,roi_topB,roi_bottomB,'automatic','automatic',0,'_B');
% Only for first 5 seconds (approx, to frame 125):
% saveVideoManyTracks('533',folder1,list1,roi_topB,roi_bottomB,'automatic',125,1,'_B_toframe125');
% ----------
% For both cells at same time
% roi_topC = [94 448 34 95];  roi_bottomC = [93 447 293 354];
% ----------------------
% ----------------------
% For data set cybD-mCherry &ATPase-GFp, image 498:
% folder2 = 'C:\Isabel\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\cybD-mCherry-ATPase-GFp_498';
% folder2 = 'Z:\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\cybD-mCherry-ATPase-GFp_498';
% list2 = [1 2 3 4 5 7 8 9 10 11 12 14 15 16 18 19 22 23 25 26 27 28 29 30 31 32 33 35 36 37 38 39 40 41 44 45 47 48 50 51 52 53 54 55 56 58 59 61 62 63];
% saveVideoManyTracks('498',folder2,list2,'full_image','full_image','automatic','automatic',0,'');
% roi498_top = [325 367 52 108];
% roi498_bottom = [323 365 312 368];
% saveVideoManyTracks('498',folder2,list2,roi498_top,roi498_bottom,'automatic',150,1,''); % forcing last frame to be frame 150.
% From 9seconds (frame 225) to 14 seconds (frame 350), forcing start and end:
% saveVideoManyTracks('498',folder2,list2,roi498_top,roi498_bottom,225,350,1,'');
% saveVideoManyTracks('498',folder2,list2,roi498_top,roi498_bottom,'automatic',100,1,'upToFrame100');
% ----------------------
% ----------------------
% For data set cybD-mCherry &ATPase-GFp, image 518:
% folder3 = 'Z:\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\cybD-mCherry-ATPase-GFp_518';
% list3 = [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22 23 24 25 26 29 30 31 32 36 37 38 40];
% roi518a_top = [228 283 16 74];
% roi518a_bottom = [228 283 256+16 256+74];
% roi518b_top = [181 270 172 235];
% roi518b_bottom = [181 270 256+172 256+235];
% saveVideoManyTracks('518',folder3,list3,roi518a_top,roi518a_bottom,150,1,'_a');
% saveVideoManyTracks('518',folder3,list3,roi518b_top,roi518b_bottom,150,1,'_b');
% -----------------------
% ----------------------
% For data set 8a: cydB-mCherry & GFPuv4-nuoF, image 484:
% folder4 = 'Z:\ExperimData\HeikoData\cydB-mCherry & GFPuv4-nuoF\cydB-mCherry & GFPuv4-nuoF TIRF\473nm770uW-561nm280uW\good\cydB-mCherry-GFPuv4-nuoF_484';
% list4 = [1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 34 35 36 37 38 39 40 41 42];
% roi484_top = [250 310 40 117];
% roi484_bottom = [250 310 300 377];
% saveVideoManyTracks('484',folder4,list4,roi484_top,roi484_bottom,100,1,'');



%% PARAMETERS

% Frame rate to show video and for the saved .avi video file:
framesPerSecond = 5;
% 25 frames per second corresponds to the experimental 40ms between frames
% (Oxphos data, Isabel). Use 5 for tracks to be more clearly visible.


%% Regions of interest:

if ~strcmp(roi_top,'full_image') && ~strcmp(roi_bottom,'full_image')
    % Top channel roi:
    xleft_top = roi_top(1);
    xright_top = roi_top(2);
    ytop_top = roi_top(3);
    ybottom_top = roi_top(4);
    % Bottom channel roi:
    xleft_bottom = roi_bottom(1);
    xright_bottom = roi_bottom(2);
    ytop_bottom = roi_bottom(3);
    ybottom_bottom = roi_bottom(4);
    
    % Error control:
    % Make sure that regions of interest for top and bottom channels are of the
    % same size:
    if (xright_top-xleft_top)~=(xright_bottom-xleft_bottom) || (ybottom_top-ytop_top)~=(ybottom_bottom-ytop_bottom)
        error('Make sure that input regions of interest for top and bottom channels are of the same size!')
    end
end


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.


%% Create output directory for saving video (.avi) results:

initial_folder_path = cd; % current initial directory.

if saveAvi == 1
    % Make new folder (new directory) to save trajectory analysis results:
    output_folder_name = strcat('manyTracksVideoAvi');
    warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
    % Create new folder for outputs inside current folder:
    mkdir(output_folder_name); % make new directory.
    output_folder_path = strcat(initial_folder_path,'\',output_folder_name);
end


%% Read trajectory data (excel files):

cd(folderWithTrajs); % move into directory containing tracks to plot.

% Loop through tracks to plot (traj numbers in input list_trajNumbers):
for i=1:length(list_trajNumbers)
    
    % Example of traj excel filename:
    % 'GFP-nuoF-mCherry-sdhC_533_traj1.xls'.
    
    % Find paths in current folder which contain 'image_label' string, the traj number, and the string '.xls':
    trajXlsPath0 = dir(strcat('*',image_label,'_traj',num2str(list_trajNumbers(i)),'.xls')); % Trajectory data (there is one excel file per trajectory in the input folder folderWithTrajs).
    % Error control:
    if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
        error('Check you are in the correct directory and run again. No .xls file found for that image number and trajectory number. Make sure image number is in between quotes ''.');
    end    
    if length(trajXlsPath0)>1
        disp(['Error for trajectory number ',num2str(list_trajNumbers(i)),' ->'])
        error('More than one excel file found for that trajectory number.')
    end
    results(i).trajXlsPath = trajXlsPath0.name; % store path of trajectory in structure "results", which has as many elements as tracks.
    
    % Import the data in the sheet named 'Track info':
    [numeric0,txt0,raw0] = xlsread(trajXlsPath0.name,'Track info');
    % Turn imported data from excel file into a structure where parameter names are fieldnames in the structure:
    str_TrackInfo = cell2struct(raw0(:,2),raw0(:,1),1);
    results(i).FirstTrajFrame = str_TrackInfo.FirstTrajFrame; 
    results(i).LastTrajFrame = str_TrackInfo.LastTrajFrame;
    results(i).NumDataPoints = str_TrackInfo.NumDataPoints;
    results(i).TopOrBottom = str_TrackInfo.TopOrBottom; % 'top' or 'bottom' region of image (colour channel, red or green);
    
    % Import the data in the sheet named 'Track data':
    [numeric,txt,raw] = xlsread(trajXlsPath0.name,'Track data');
    % Turn imported data from excel file into a structure:
    str_TrackData = struct; % create empty structure to fill up.
    for j = 1:length(txt)
        str_TrackData = setfield(str_TrackData, txt{j}, numeric(:,j)); % create field in the structure.
        % Each field in the structure contains a column vector with the data.
    end
    results(i).frame = str_TrackData.frame; % column vector with frame number for that track
    results(i).xvalues0 = str_TrackData.xvalues0; % column vector with original x-values for track.
    results(i).yvalues0 = str_TrackData.yvalues0; % column vector with original y-values for track.
end

cd(initial_folder_path); % return to initial directory.

% The result of the previous loop is the structure "results", which has as
% many elements as trajectories to plot (length(list_trajNumbers)), each
% with fields: trajXlsPath, FirstTrajFrame, LastTrajFrame, NumDataPoints,
% TopOrBottom, frame, xvalues0, yvalues0.

% First and last frames to show on video with all tracks:
if strcmp(forced_end_frame,'automatic')
    last_frame = max([results.LastTrajFrame]); % max of list of last frames of all tracks.
else
    last_frame = forced_end_frame; % force end of video to this particular frame number.
end

if strcmp(forced_start_frame,'automatic')
    first_frame = min([results.FirstTrajFrame]); % min of list of first frames of all tracks.
else
    first_frame = forced_start_frame; % force start of video to this particular frame number.
end


%% Show video and save it (if input "saveAvi" is equal to 1):

clear mex % close all open .avi files (this avoids future errors).

if saveAvi == 1
    % Move into output folder to save videos:
    cd(output_folder_path);
end

fig = figure('Tag','Data video','units','inches','position',[12 4 6 6]); % Figure number 2.
% 'position' vector is [left, bottom, width, height].
% left, bottom control the position at which the window appears when it pops.
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')

if saveAvi == 1
    video_filename = strcat('ManyTracksVideo_',image_label,output_label,'.avi'); % filename of video.
    % Create and open avi file for saving frames onto it later:    
    mov = avifile(video_filename,'fps',framesPerSecond,'compression','None'); % this saves the .avi file in the current directory automatically. Movie with circle overlaid onto detected spot.
    % 'fps' is the frame rate for the saved video, 25 frames per second
    % corresponds to the experimental 40ms between frames. 
    % 'compression' codecs can be 'Indeo3', 'Indeo5'(default), 'Cinepak', 'MSVC',
    % 'RLE' or 'None', only 'None' works.    
end


% Loop through frames:
% Start at the very first frame in all trajectories and end at the very
% last frame in all trajs.

for n = first_frame:last_frame
  
    if strcmp(roi_top,'full_image') && strcmp(roi_bottom,'full_image') % for full image.
        frame = image_data(n).frame_data; % extract frame data, stored in the field 'frame_data'.        
    else % with regions of interest:       
        frame0 = image_data(n).frame_data; % extract frame data.
        frame_top = mat2gray(frame0(ytop_top:ybottom_top,xleft_top:xright_top));
        frame_bottom = mat2gray(frame0(ytop_bottom:ybottom_bottom,xleft_bottom:xright_bottom));
        frame = [frame_top; frame_bottom]; % combine both images.
    end
    
    % Show frame:
    imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).   
    hold on;
    % Plot trajectories overlaid on frame if there are any for this frame number:
    for q = 1:length(list_trajNumbers) % Loop through list of trajectories to plot.
         
        % See if current frame, "n" is part of that trajectory:
        list_frames_inTrack = results(q).frame; % column vector, list of frame numbers in track "q".               
        logic_vector = ismember(list_frames_inTrack,n); % vector with 1 at position of "n" in vector "list_frames_inTrack".
        
        if sum(logic_vector)~=0 % if the current frame number, "n", is part of that trajectory:
            
             xvector = results(q).xvalues0; % column vector.
             yvector = results(q).yvalues0; % column vector.
             % Extract position of bright spot for frame "n" in trajectory
             % given by list_trajNumbers(q):
             xposition0 = xvector(logic_vector==1);
             yposition0 = yvector(logic_vector==1);
            
            % Decide color for plotting circle around bright spot:
            if strcmp(results(q).TopOrBottom,'top')
                colour_forCircle = 'r';
            elseif strcmp(results(q).TopOrBottom,'bottom')
                colour_forCircle = 'g';
            end
            % Plot circle around spot:
            if strcmp(roi_top,'full_image') && strcmp(roi_bottom,'full_image') % for full image.
                new_x_position = xposition0;
                new_y_position = yposition0;
            else % with regions of interest:
                if strcmp(results(q).TopOrBottom,'top')
                    new_x_position = xposition0-xleft_top+1;
                    new_y_position = yposition0-ytop_top+1;
                elseif strcmp(results(q).TopOrBottom,'bottom')
                    new_x_position = xposition0-xleft_bottom+1;
                    new_y_position = yposition0-ytop_bottom+1+(ybottom_top-ytop_top+1);
                end
            end
            
            plot(new_x_position,new_y_position,'o','Color',colour_forCircle,'MarkerSize',10)
        end
                
    end
    
    hold off;
    % pause(1/framesPerSecond) % 40ms pause in between shown images (frame rate of saved video is determined by avifile parameter 'fps').
    pause(0.2)
    
    if saveAvi == 1
        % Save frame with circle overlaid, to the other video, mov:
        F = getframe(gca);
        mov = addframe(mov,F);        
    end
            
    close(findobj('Tag','Data video')); % close video figure;
    
end


if saveAvi == 1
    % Close avi file:
    mov = close(mov);
    % Return to previous folder:
    cd(initial_folder_path);
end