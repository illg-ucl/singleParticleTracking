function results = fullCellIntensity(dataSet_label,image_label,start_frame,tsamp,roi,output_label) 
%
% Created by Isabel Llorente-Garcia, June 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% 
% This function looks at the integrated intensity of a whole bacteria cell
% (subtracting background) and its decay over time to try and possibly find
% the photobleaching time constant for a given fluorophore.
% Note that to run this function the current directory needs to be that
% which contains the video (image sequence) file identified by input
% image_label.
%
% Inputs: 
% - dataSet_label: string which labels the data set. Eg. 'ATPase-GFP_'. This
% will form part of the name of the folder created to save results for a given image.
% - image_label: string which identifies the image sequence, eg. '110' (string in video name before '.sif').
% - start_frame: number corresponding to first frame to go through.
% - 'tsamp' is the sampling time, ie., the time between frames in seconds.
%    Eg. for 40ms per frame, tsamp is 0.04. 
% - roi: initial region of interest: enter as vector [xleft xright ytop ybottom] (in pixels). 
% Eg. [1 512 256 512] for the lower half of the image,
% Eg. [1 512 1 256] for the top half of the image,
% Eg. [1 512 1 512] for the full image.
% - output_label: string which can be '_cell1', '_cell2', 'bis', ' ', etc... to
% add at end of name of output excel file. 
%
% Output: 
% The output, "results", is a cell array with two elements.
% The first element is a structure containing the total cell intensity
% trace (time vector relative to start of image sequence and total cell
% intensity after bgnd subtraction).
% The second element in the cell array is a structure containing the fit
% results.
%
% Example:  
% w110 = fullCellIntensity('ATPase-GFP','110',37,0.04,[1 512 257 512],' ');


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.

%% PARAMETERS

% PARAMETER r:
% Number of frames to average over when calculating a frame average in
% order to get a Signal Mask to distinguish cell region from background
% region. Note it has to be at most equal to the total no. of frames minus
% start_frame. Sometimes it help to only use the first frames to get a
% better cell signal mask.
% CHECK: check if only a few first frames (20) needed or more:
r = numFrames-start_frame; % number of frames to average over, starting from start_frame. (20, numFrames-start_frame)

% CHECK:
% PARAMETERs for default cell region around whole cell:
xsize_default_cell_region = 140; % horizontal size of region around cell. (70-140) 
ysize_default_cell_region = 140; % horizontal size of region around cell. (70-140)

% How many pixels to take from the cell edge for a rectangle around the cell, for background:
pix_away = 8; % default is 10 pixels.


%% Initial region of interest, input roi;
xleft_0 = roi(1);
xright_0 = roi(2);
ytop_0 = roi(3);
ybottom_0 = roi(4);


%% Create new directory for saving trajectory-analysis results for this image sequence:

% Make new folder (new directory) to save trajectory analysis results:
new_folder_name = strcat(dataSet_label,'_wholeCell_I');
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(new_folder_name); % make new directory.

% Name of output excel file with results:
output_filename = strcat(dataSet_label,'_',image_label,output_label,'.xls'); % name of excel file to save to.

%% Calculate frame average:

% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul = zeros(frame_Ysize,frame_Xsize);

for k = start_frame:numFrames  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul = frame_accumul + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg = frame_accumul/(numFrames-start_frame+1); % values are average ones, real ones (not between 0 and 1).
frame_avg_roi = frame_avg(ytop_0:ybottom_0,xleft_0:xright_0);


%% Calculate Signal Mask to distinguish cell region from background region:
% Use frame average (of all frames or first frames only) to calculate signal mask. 

% Calculate average of first "r" frames:
% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul_0 = zeros(frame_Ysize,frame_Xsize);
% r is the number of frames to average over, starting from start_frame.
% See PARAMETERS section.

for k = start_frame:start_frame+r  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul_0 = frame_accumul_0 + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg_0 = frame_accumul_0/(r+1);
frame_avg_0_roi = frame_avg_0(ytop_0:ybottom_0,xleft_0:xright_0);
frame_avg_0_roi_Gray = mat2gray(frame_avg_0_roi); % The input to function "getCellMaskAndBoundary" needs to be a grayscale image:

% Threshold the selected roi of the full image to find a preliminary signal mask. 
% A local threshold value within a selected region given by the second
% input parameter (full image here) in getCellMaskAndBoundary2 is used.
% Get SignalMask to know where cells are, to distinguish cells from background:
[SignalMask_0 CellBoundaryMask_0] = getCellMaskAndBoundary2(frame_avg_0_roi_Gray,[1 size(frame_avg_0_roi_Gray,2) 1 size(frame_avg_0_roi_Gray,1)]); 
% % SignalMask is a matrix with 1 at positions where cells are and 0 at
% % background.

% Centre of mass:
% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos. For the whole image:
[Xpos,Ypos] = meshgrid(1:size(SignalMask_0,2),1:size(SignalMask_0,1));
Cx = SignalMask_0.*Xpos; % Positions of ones within signal mask.
Cy = SignalMask_0.*Ypos; % Positions of ones within signal mask.
Cx_1 = Cx(Cx~=0); % extract non-zero values. X values of cell region for selected region of interest.
Cy_1 = Cy(Cy~=0); % extract non-zero values. Y values of cell region for selected region of interest.
% Centre of mass of preliminary cell region:
centre_of_mass_X_0 = mean(Cx_1);
centre_of_mass_Y_0 = mean(Cy_1);

% Default cell region: it is a square of size (xsize_default_cell_region,
% ysize_default_cell_region) around the centre of mass found:
xleft = round(centre_of_mass_X_0 - xsize_default_cell_region/2);
xright = round(centre_of_mass_X_0 + xsize_default_cell_region/2);
ytop = round(centre_of_mass_Y_0 - ysize_default_cell_region/2);
ybottom = round(centre_of_mass_Y_0 + ysize_default_cell_region/2);

% Error control:
% if box for default_cell_region too close to image edge, limit value to image size:
xleft = max(1,xleft); % this makes xleft >=1.
xright = min(xright,size(SignalMask_0,2)); % this makes xright<=size(SignalMask_0,2). 
ytop = max(1,ytop); % this makes ytop >=1.
ybottom = min(ybottom,size(SignalMask_0,1)); % this makes ybottom<=size(SignalMask_0,1). 

% Default cell region corners and centre (for calculation of local cell coordinates) for plot:
def_cell_region_x = [xleft,xleft,xright,xright]; % corners of default cell region, x values.
def_cell_region_y = [ytop,ybottom,ytop,ybottom]; % corners of default cell region, y values.

% Plot frame average, preliminary cell signal mask and default cell region:
% figure(200) % create figure number 200.
h = figure('Tag','cell_cords','position',[100 300 1500 800]); 
% 'position' vector in figure is [left, bottom, width, height].

subplot(2,2,1);
imshow(frame_avg,[],'Border','tight'); % display frame average (full image) scaled between its min and max values ([]).

subplot(2,2,3);
imshow(SignalMask_0,[],'Border','tight'); % display preliminary cell signal mask around cell region in ROI.
title('Preliminary cell mask');
hold on;
% Plot centre of mass of default cell region in red:
plot(centre_of_mass_X_0,centre_of_mass_Y_0,'+r','MarkerSize',3) 
% Plot default cell region corners in green:
plot(def_cell_region_x,def_cell_region_y,'+g','MarkerSize',3) 
hold off;


% CHECK: give option of user input or not. Comment or uncomment as you wish.
% use_default_cell_region = input('Use plotted cell region for local cell coordinates (1 for yes, 0 for no)? :');
use_default_cell_region = 1;
local_cell_region = [xleft xright ytop ybottom];
% CHECK: comment out some of the following or not:
% Particular code for a given image with two cells, one left, one right:
if use_default_cell_region == 0
    % Uncomment one of the following options:
    % local_cell_region = input('Enter region around cell of interest as
    % [xleft xright ytop ybottom] (for local cell coordinates):'); % request user input
    % --------------------------------
    % local_cell_region = [205 264 152 209];
    % --------------------------------
%     if centre_of_mass_X_0 < 250 % cell on the left
%         local_cell_region = [192 244 149 200];
%     else % cell on the right
%         local_cell_region = [289 349 17 56];
%     end
    % --------------------------------
    if centre_of_mass_X_0 < 68
        local_cell_region = [21 73 62 105];
    else % top cells
        if centre_of_mass_X_0 < 150
            local_cell_region = [70 113 116 149];
        else
            local_cell_region = [200 251 139 202];
        end
    end
    % --------------------------------
    xleft = local_cell_region(1);
    xright = local_cell_region(2);
    ytop = local_cell_region(3);
    ybottom = local_cell_region(4);
end


% CHECK:
% Use a local threshold within the selected cell region to threshold the
% image and find the cell-signal mask: use getCellMaskAndBoundary2.m instead
% of getCellMaskAndBoundary.m.
[SignalMask CellBoundaryMask] = getCellMaskAndBoundary2(frame_avg_0_roi_Gray,local_cell_region); 
% frame_avg_0_Gray is the complete greyscale image.
BgndMask = ~SignalMask; % background mask.

% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos. For the whole image:
[Xpos,Ypos] = meshgrid(1:size(SignalMask,2),1:size(SignalMask,1));
Cx = SignalMask.*Xpos; % Positions of ones within signal mask.
Cy = SignalMask.*Ypos; % Positions of ones within signal mask.
Cx_1 = Cx(Cx~=0); % extract non-zero values. X values of cell region for selected region of interest.
Cy_1 = Cy(Cy~=0); % extract non-zero values. Y values of cell region for selected region of interest.
% Centre of mass of final cell region:
centre_of_mass_X = mean(Cx_1);
centre_of_mass_Y = mean(Cy_1);
minX_cellMask = min(Cx_1);
maxX_cellMask = max(Cx_1);
minY_cellMask = min(Cy_1);
maxY_cellMask = max(Cy_1);

% Small box region around cell, delimiting edges of box:
ytop_box = minY_cellMask - pix_away;
ybottom_box = maxY_cellMask + pix_away;
xleft_box = minX_cellMask - pix_away;
xright_box = maxX_cellMask + pix_away;

box_x = [xleft_box,xleft_box,xright_box,xright_box]; % corners of box, x values.
box_y = [ytop_box,ybottom_box,ytop_box,ybottom_box]; % corners of box, y values.

% Error control:
% if box for default_cell_region too close to image edge, limit value to image size:
xleft_box = max(1,xleft_box); % this makes xleft >=1.
xright_box = min(xright_box,size(SignalMask,2)); % this makes xright<=size(SignalMask,2). 
ytop_box = max(1,ytop_box); % this makes ytop >=1.
ybottom_box = min(ybottom_box,size(SignalMask,1)); % this makes ybottom<=size(SignalMask,1). 
% Use a small rectangular box around the cell-region, "pix_away" pixels
% away from the cell edge:
SignalMask_box = SignalMask(ytop_box:ybottom_box,xleft_box:xright_box);
BgndMask_box = BgndMask(ytop_box:ybottom_box,xleft_box:xright_box);
% pos_bgnd_box = find(BgndMask_box==1); % positions of ones in bgnd mask. (used later)

% Plot final cell signal mask after using a local threshold (substituting previous plot):
subplot(2,2,3);
imshow(SignalMask,[],'Border','tight'); % display preliminary cell signal mask around cell region.
title('Final cell mask (uses local threshold)');

subplot(2,2,2);
imshow(frame_avg_roi,[],'Border','tight'); % display preliminary cell signal mask around cell region.
title('Frame average, selected roi around cell');
hold on;
% Plot centre of mass of default cell region in red:
plot(centre_of_mass_X,centre_of_mass_Y,'+r','MarkerSize',3) 
% Plot box corners in green:
plot(box_x,box_y,'+g','MarkerSize',3) 
hold off;

% % To check graphically on frame average:
% subplot(2,2,4);
% cell_bgnd = BgndMask.*frame_avg_roi;
% % imshow([cell_signal; cell_bgnd],[],'Border','tight'); 
% imshow(cell_bgnd,[],'Border','tight'); 
% title('Cell signal and background');

% % Enter coordinates [x y] of centre of region to use as background mask in
% % each frame:
%  centre_point_regionForBgnd = input('Enter coordinates [x y] for centre of background region away from cell:');
% % Translate signal mask to a side to find out background:
% se = translate(strel(1), [0 1.2*xsize_default_cell_region]);
% displaced_signal_mask = imdilate(SignalMask,se);
% % imshow(displaced_signal_mask);


%% Loop through each frame in image sequence:

Itotal_cell = zeros(numFrames-start_frame+1,1); % initialise column vector of integrated intensities for each frame.
n = 1; % start counter for used frames.

for k = start_frame:numFrames  % loop through frames.
    
    % Get frame data:
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    
    % Take roi (input parameter) of frame:
    frame_roi = frame(ytop_0:ybottom_0,xleft_0:xright_0);
    % Use a small rectangular box around the cell-region, "pix_away" pixels
    % away from the cell edge:
    frame_roi_box = frame_roi(ytop_box:ybottom_box,xleft_box:xright_box);

%     % Use cell and bgnd masks:
%     cell_signal_box = SignalMask_box.*frame_roi_box;
%     cell_bgnd_box = BgndMask_box.*frame_roi_box;
    
    % Total integrated intensities within the little box:
    % Itotal_bgnd = sum(frame_roi_box(BgndMask_box==1)); % total background intensity.
    Iavg_bgnd = median(frame_roi_box(BgndMask_box==1)); % we use the median as the average bgnd intensity value per pixel.
    
    % Total integrated intensity in cell region before background subtraction:    
    Itotal_cell_0 = sum(frame_roi_box(SignalMask_box==1)); 
    % Subtract background:
    numPixInCellMask = length(frame_roi_box(SignalMask_box==1)); % no. of pixels in cell region
    % Store intensity value in column vector:
    Itotal_cell(n) = Itotal_cell_0 - numPixInCellMask*Iavg_bgnd; % Total cell intensity after bgnd subtraction.
    
    n = n+1;
   
end

frame_nums = (start_frame:numFrames)'; % column vector with frame numbers
time_vector = tsamp.*(frame_nums-start_frame); % time relative to start of image sequence (start_frame).

% Store in result structure:
results_Itrace.time_rel = time_vector;
results_Itrace.Itotal_cell = Itotal_cell;


%% %% Fit intensity data to an exponentially decaying function (with no offset):

tforFit = time_vector; % time in seconds relative to start of image sequence.
IforFit = Itotal_cell; % total cell intensity after background subtraction.

exp_no_offset = fittype('I0*exp(-t/tau)','independent','t'); % define exponential funtion to fit to, with 't' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = max(IforFit); % since I has been normalised.
guess_tau = 2; % in seconds;
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 0]; % Lower bounds for fit parameters. In order: I0, tau.
options.Upper = [Inf 100]; % Upper bounds for fit parameters. In order: I0, tau.

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
disp('Intensity vs time exponential fit (no offset) result: ') 
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
guess_I0 = max(IforFit); 
guess_Ioffset = 0;
guess_tau = 2; 
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_Ioffset guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 -0.2*max(IforFit) 0]; % Lower bounds for fit parameters. In order: I0, Ioffset, tau.
options.Upper = [Inf 0.2*max(IforFit) 100]; % Upper bounds for fit parameters. In order: I0, Ioffset, tau.

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
results_I_fits.numPixInCellMask = numPixInCellMask; % no. of pixels in cell region


%% Plot total cell intensity trace and its exponential fits:
subplot(2,2,4);
plot(time_vector,Itotal_cell,'.k')
hold on;
plot(fit_result_I,'b'); % plot exponential fit (no offset) as BLUE line.
plot(fit_result_I_wo,'r'); % plot exponential fit (with offset) as RED line.
legend('data','exp fit-no offset','exp fit-with offset');
% legend('hide');
xlabel('time from start (s)');
ylabel('IwholeCell');
title('Total cell intensity');
hold off;


%% Save figure as png in folder created

figName = strcat(dataSet_label,'_',image_label,output_label); % name of figure file to save to. 
saveFigurePNG(new_folder_name,figName)

%% Output results:
% The output, "results", is a cell array with two elements.
% The first element is a structure containing the total cell intensity
% trace (time vector relative to start of image sequence and total cell
% intensity after bgnd subtraction).
% The second element in the cell array is a structure containing the fit
% results.
results{1} = results_Itrace;
results{2} = results_I_fits;


%% Save excel file with results:

% Move into folder previously created to save traj analysis results:
cd(new_folder_name); % move into that directory.

% Save results in an excel file in the same folder where the
% analysis plot .png file has been saved, new_folder_name:
dataForSheet1 = [fieldnames(results{1})'; num2cell(cell2mat(struct2cell(results{1})'))]; % vectors
dataForSheet2 = [fieldnames(results{2}) struct2cell(results{2})]; % numbers
warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.

xlswrite(output_filename,dataForSheet1,'Itrace'); % write data to sheet 'Itrace' in excel file.
xlswrite(output_filename,dataForSheet2,'Results I fits'); % write data to sheet 'Results I fits' in excel file.

cd('..'); % go back to previous directory.

