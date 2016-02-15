% Choose one of the following options: 
% image_label = '498'; % for image 498.
image_label = '518'; % for image 518.

% Create new folder for outputs inside current folder:
output_folder_name = 'singleFrames';
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(output_folder_name); % make new directory.

if strcmp(image_label,'498')
    % For the cell on the left:
    start_frame = 52; % frames to save as individual images (roi)
    end_frame = 59; % frames to save as individual images (roi)
    mid_frame = round(mean(start_frame:end_frame));
    % Region of interest (roi) around cell. Use region of interest for all frames:
    roi_top = [315 380 40 120]; % [xleft xright ytop ybottom] for cell on the left, top channel.
    roi_bottom = [310 375 300 380]; % [xleft xright ytop ybottom] for cell on the left, bottom channel.
    % For frame average (see "calculateFrameAverages.m"):
    start_frame_forAvg = 7;
    end_frame_forAvg = 500;
    
else
    if strcmp(image_label,'518')
        % For the cell on the right:
        start_frame = 8; % frames to save as individual images (roi)
        end_frame = 15; % frames to save as individual images (roi)
        mid_frame = round(mean(start_frame:end_frame));
        % Region of interest (roi) around cell. Use region of interest for all frames:
        roi_top = [230 280 15 75]; % [xleft xright ytop ybottom] for cell on the left, top channel.
        roi_bottom = [225 275 273 333]; % [xleft xright ytop ybottom] for cell on the left, bottom channel.
        % For frame average (see "calculateFrameAverages.m"):
        start_frame_forAvg = 7;
        end_frame_forAvg = 500;
    end
end

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);

% Use max and min of mid_frame to plot all frames with same contrast:
mid_frame_data = image_data(mid_frame).frame_data; % extract matrix data for first frame.
mid_frame_data = double(mid_frame_data);
mid_frame_data_top = mid_frame_data(roi_top(3):roi_top(4),roi_top(1):roi_top(2)); % select top roi.
mid_frame_data_bottom = mid_frame_data(roi_bottom(3):roi_bottom(4),roi_bottom(1):roi_bottom(2)); % select bottom roi.
% Get max and min of top and bottom regions of interest:
mid_frame_top_min = min(min(mid_frame_data_top)); % min top channel.
mid_frame_top_max = max(max(mid_frame_data_top)); % max top channel.
mid_frame_bottom_min = min(min(mid_frame_data_bottom)); % min bottom channel.
mid_frame_bottom_max = max(max(mid_frame_data_bottom)); % max bottom channel.

for i = start_frame:end_frame
    
    % Read image-sequence file:
    frame = image_data(i).frame_data; % extract matrix data for first frame.
    frame = double(frame);
    frame_top = frame(roi_top(3):roi_top(4),roi_top(1):roi_top(2));
    frame_bottom = frame(roi_bottom(3):roi_bottom(4),roi_bottom(1):roi_bottom(2));
    disp(['frame number: ',num2str(i)]) % print frame number to Command Window.
    
    % Show top channel roi:
    figure;
    imshow(frame_top,[mid_frame_top_min mid_frame_top_max])
    % Save current figure as .png in output folder and close image:
    figName_top = strcat('Image',image_label,'frame',num2str(i),'_','top'); % name of figure file to save to.
    saveFigurePNG(output_folder_name,figName_top)
    
    % Show bottom channel roi:
    figure;
    imshow(frame_bottom,[mid_frame_bottom_min mid_frame_bottom_max])
    % Save current figure as .png in output folder and close image:
    figName_bottom = strcat('Image',image_label,'frame',num2str(i),'_','bottom'); % name of figure file to save to.
    saveFigurePNG(output_folder_name,figName_bottom)
    
end

%% Now frame averages in same roi:

clear('image_data'); % clear variable to avoid memory problems.
A = frameAverage(image_label,start_frame_forAvg,end_frame_forAvg,0,0); % calculate frame average (do not plot, do not save (0,0)).

A2_top = A(roi_top(3):roi_top(4),roi_top(1):roi_top(2)); % select top roi.
A2_bottom = A(roi_bottom(3):roi_bottom(4),roi_bottom(1):roi_bottom(2)); % select bottom roi.

% Show top channel roi:
figure;
imshow(A2_top,[])
% Save current figure as .png in output folder and close image:
figName_top2 = strcat('Image',image_label,'frameAvgROI_','top'); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName_top2)

% Show bottom channel roi:
figure;
imshow(A2_bottom,[])
% Save current figure as .png in output folder and close image:
figName_bottom2 = strcat('Image',image_label,'frameAvgROI_','bottom'); % name of figure file to save to.
saveFigurePNG(output_folder_name,figName_bottom2)
