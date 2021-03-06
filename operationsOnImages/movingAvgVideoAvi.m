function movingAvgVideoAvi(image_label,start_frame,end_frame,numFramesToAverageOver,roi_top,roi_bottom,scaling_option,saveAvi,output_label)
% 
% Created by Isabel Llorente-Garcia, August 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Create a moving average video of a given image sequence and save (or not)
% as an .avi video file.
%
% Inputs:
% - image_label: label for image sequence of choice, string such as '498'.
% - start_frame and end_frame: do the moving average from start_frame to
% end_frame. Note that "end_frame" can be the string 'end', ie., the last
% frame.
% - numFramesToAverageOver: width of window to average over. It must be an
% odd number, since we average a given frame with a number of frames [(numFramesToAverageOver-1)/2] before
% and after it.
% - roi_top and roi_bottom: "regions of interest" around the cell, for top
% and bottom channels, regions to show on video. "roi" is a row vector
% [xleft xright ytop ybottom]. If equal to 'full_image', then the whole image is used.
% - scaling_option: option (1, 2 or 3) for re-scaling images. Option 1:
% re-scaling takes place after top and bottom channels
% (original values, unmodified) are put together; so top and bottom
% channels are scaled together and plotted between the min and max in the whole image, frame by
% frame. Option 2: each channels is first re-scaled between its mean and
% max independently of the other channel, and they are then put together.
% Option 3: use a given frame number to find the min and max to use for
% re-scaling each channel independently, and all frames to those same min
% and max values.
% - saveAvi: save as .avi videofile (1) or do not save (0).
% - output_label: string for labeling the output .avi file, for different
% cells or roi's in one image, for instance.
%
% Example of how to run this function:
% movingAvgVideoAvi('498',7,'end',3,'full_image','full_image',2,1,'A') to save avi video, for
% full image.
% movingAvgVideoAvi('498',7,'end',3,[300 390 40 130],[300 390 256+40 256+130],2,0,'') not saving .avi, for a roi.
% ---------------
% roi498_top2 = [325 367 52 108]; roi498_bottom2 = [323 365 312 368];
% movingAvgVideoAvi('498',20,270,3,roi498_top2,roi498_bottom2,3,0,'')


%% PARAMETERS

% Frame rate to show video and for the saved .avi video file:
framesPerSecond = 25;
% 25 frames per second corresponds to the experimental 40ms between frames.

% Frame to use for scaling option 3: re-scaling all frames between the min
% and max in that frame number, and for top and bottom channels
% independently:
% frame_num_for_scaling = (start_frame+end_frame)/2;
% frame_num_for_scaling = start_frame+round((end_frame-start_frame)/2);
frame_num_for_scaling = start_frame+round((end_frame-start_frame)/2.5);


%% Error control: numFramesToAverageOver must be an odd number
if mod(numFramesToAverageOver,2)==0 % if numFramesToAverageOver is an even number
    return % exit program.
end

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
end

%% Create output directory for saving video (.avi) results:

if saveAvi == 1
    % Make new folder (new directory) to save trajectory analysis results:
    initial_folder_path = cd;
    % Create new folder for outputs inside current folder:
    output_folder_name = strcat('movingAvgVideosAvi');
    warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
    mkdir(output_folder_name); % make new directory.
    output_folder_path = strcat(initial_folder_path,'\',output_folder_name);
end


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.

if strcmp(end_frame,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   end_frame = numFrames;
end

clear mex % close all open .avi files (this avoids future errors).

if saveAvi == 1
    % Move into output folder to save videos:
    cd(output_folder_path);
end


fig = figure('Tag','Data video','units','inches','position',[12 4 6 6]); % Figure number 2.
% 'position' vector is [left, bottom, width, height].
% left, bottom control the position at which the window appears when it
% pops.
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],'NextPlot','replace','Visible','off')

if saveAvi == 1
    video_filename = strcat('MovingAvgVideo_',image_label,output_label,'.avi'); % filename of video.
    % Create and open avi file for saving frames onto it later:    
    mov = avifile(video_filename,'fps',framesPerSecond,'compression','None'); % this saves the .avi file in the current directory automatically. Movie with circle overlaid onto detected spot.
    % 'fps' is the frame rate for the saved video, 25 frames per second
    % corresponds to the experimental 40ms between frames. 
    % 'compression' codecs can be 'Indeo3', 'Indeo5'(default), 'Cinepak', 'MSVC',
    % 'RLE' or 'None', only 'None' works.    
end

% SCALING OPTION 3:
% For scaling option 3, find out min and max for re-scaling:
m = frame_num_for_scaling; % frame number used for re-scaling all frames, top and bottom independently.
frame0 = image_data(1).frame_data; % (used to get size.)
if scaling_option == 3
    % Find out rolling average within average window around frame "m":
    if strcmp(roi_top,'full_image') && strcmp(roi_bottom,'full_image') % for full image.
        frame_accumul = single(zeros(size(frame0))); % initialise empty matrix of same size as each image frame.
        for k = (m-(numFramesToAverageOver-1)/2):(m+(numFramesToAverageOver-1)/2)
            % Get single frame data:
            frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
            frame_accumul = frame_accumul + frame;
        end
        frame_to_plot = frame_accumul/numFramesToAverageOver;
        % frame_accumul/numFramesToAverageOver is the frame average of
        % numFramesToAverageOver frames around frame m.
        scale_min = min(min(frame_to_plot));
        scale_max = max(max(frame_to_plot));
    else % with regions of interest:
        frame_accumul_top = single(zeros(ybottom_top-ytop_top+1,xright_top-xleft_top+1)); % initialise empty matrix of same size as each roi.
        frame_accumul_bottom = single(zeros(ybottom_bottom-ytop_bottom+1,xright_bottom-xleft_bottom+1)); % initialise empty matrix of same size as each roi.
        for k = (m-(numFramesToAverageOver-1)/2):(m+(numFramesToAverageOver-1)/2)
            % Get single frame data:
            frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.  
            % Separate top and bottom regions of interest:
            frame_top = frame(ytop_top:ybottom_top,xleft_top:xright_top);
            frame_bottom = frame(ytop_bottom:ybottom_bottom,xleft_bottom:xright_bottom);    
            % Add up and accumulate frames:
            frame_accumul_top = frame_accumul_top + frame_top;
            frame_accumul_bottom = frame_accumul_bottom + frame_bottom;           
        end
        frame_to_plot_top = frame_accumul_top/numFramesToAverageOver;
        frame_to_plot_bottom = frame_accumul_bottom/numFramesToAverageOver;
        % find  min and max for top and bottom channels independently:
        scale_min_top = min(min(frame_to_plot_top));
        scale_max_top = max(max(frame_to_plot_top));
        scale_min_bottom = min(min(frame_to_plot_bottom));
        scale_max_bottom = max(max(frame_to_plot_bottom));        
    end
end
    


% Loop through frames:
% Start at the first frame for which I can take the whole averaging window
% and go until the last frame for which I can do the same:
for n = (start_frame+(numFramesToAverageOver-1)/2):(end_frame-(numFramesToAverageOver-1)/2)
  
    % Within the averaging window:

    if strcmp(roi_top,'full_image') && strcmp(roi_bottom,'full_image') % for full image.
        frame_accumul = single(zeros(size(frame0))); % initialise empty matrix of same size as each image frame.
        
        for k = (n-(numFramesToAverageOver-1)/2):(n+(numFramesToAverageOver-1)/2)
            % Get single frame data:
            frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
%             frame = double(frame);
            frame_accumul = frame_accumul + frame;
        end
        frame_to_plot = frame_accumul/numFramesToAverageOver;
        if scaling_option == 3
            frame_to_plot = (frame_to_plot-scale_min)/(scale_max-scale_min); % re-scale full image between min and max of frame number "frame_num_for_scaling".
        end
        
    else % with regions of interest:
        frame_accumul_top = single(zeros(ybottom_top-ytop_top+1,xright_top-xleft_top+1)); % initialise empty matrix of same size as each roi.
        frame_accumul_bottom = single(zeros(ybottom_bottom-ytop_bottom+1,xright_bottom-xleft_bottom+1)); % initialise empty matrix of same size as each roi.
        
        for k = (n-(numFramesToAverageOver-1)/2):(n+(numFramesToAverageOver-1)/2)
            % Get single frame data:
            frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
            % frame = double(frame);
            
            if scaling_option == 1 || scaling_option == 3
                frame_top = frame(ytop_top:ybottom_top,xleft_top:xright_top);
                frame_bottom = frame(ytop_bottom:ybottom_bottom,xleft_bottom:xright_bottom);
            elseif scaling_option == 2
                frame_top = mat2gray(frame(ytop_top:ybottom_top,xleft_top:xright_top));
                frame_bottom = mat2gray(frame(ytop_bottom:ybottom_bottom,xleft_bottom:xright_bottom));
            end
            % add up and accumulate frames for each channel:
            frame_accumul_top = frame_accumul_top + frame_top;
            frame_accumul_bottom = frame_accumul_bottom + frame_bottom;
        end
        % Divide by number of frames to get the average:
        frame_to_plot_top = frame_accumul_top/numFramesToAverageOver;
        frame_to_plot_bottom = frame_accumul_bottom/numFramesToAverageOver;
        if scaling_option == 3 % re-scale to the min and max of the chosen frame number:
            frame_to_plot_top = (frame_to_plot_top-scale_min_top)/(scale_max_top-scale_min_top);
            frame_to_plot_bottom = (frame_to_plot_bottom-scale_min_bottom)/(scale_max_bottom-scale_min_bottom);
        end        
        % Combine both channels together:
        frame_to_plot = [frame_to_plot_top; frame_to_plot_bottom]; 
        
    end
    
    % Show averaged frame:
    if scaling_option == 1 || scaling_option == 2
        imshow(frame_to_plot,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
    elseif scaling_option == 3        
        imshow(frame_to_plot,[0 1],'Border','tight','InitialMagnification',150); % show image scaled between 0 and 1.
    end
    
    pause(1/framesPerSecond) % 40ms pause in between shown images (frame rate of saved video is determined by avifile parameter 'fps').
    
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