function [shift_y shift_x] = overlapBF(image_label,start_frame,end_frame,exclude_y_top,exclude_y_bot,exclude_x_left,exclude_x_right,show_plots,saveplot)
% 
% Isabel Llorente-Garcia: August 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% For Bright Field Images, to find shift between top and bottom channels.
% Similar to quickGraphColocalisation.m but for BF images.
%
% Inputs:
% - image_label: such as '513', '490', etc... which corresponds to a certain .sif number of image sequence
% - start_frame and end_frame: frames to use. If end_frame is 'end', all
% frames are used.
% - exclude_y_top,exclude_y_bot, exclude_x_left,exclude_x_right: exclude
% edges of the image (black edges which appear on bright field images mess
% up the cross-correlation calculation). Exclude a number of pixels on each
% edge equal to exclude_y_top,exclude_y_bot,
% exclude_x_left,exclude_x_right.
% - show_plots: 1 for yes, 0 for no plots.
% 
% Example of how to call this function:
% [y0 x0] = overlapBF('497',1,10,20,30,20,35,1)
% [y0 x0] = overlapBF('497',1,'end',20,30,20,35,1)
% Bright field images are smaller in bites, usually only 10 frames.


% Calculate frame average between start_frame and end_frame:
fr_avg = frameAverage(image_label,start_frame,end_frame,0,0); % don't show or save.

% Separate top and bottom channels (usually red and green channels) and
% exclude edges of images.
% Also convert to grayscale (0 to 1 values) images (so final result is only
% qualitative):
Btop = mat2gray(fr_avg(1+exclude_y_top:(size(fr_avg,1)/2-exclude_y_bot),1+exclude_x_left:size(fr_avg,2)-exclude_x_right)); % top channel.
Bbot = mat2gray(fr_avg((size(fr_avg,1)/2+1+exclude_y_top):(size(fr_avg,1)-exclude_y_bot),1+exclude_x_left:size(fr_avg,2)-exclude_x_right)); % bottom channel.

% --------------------------------------------------
% Spatial cross-correlation of top and bottom images to find out the
% spatial shift between red and green (top and bottom) channels:
size_cor = 30; % how far away to go on displacement in each direction to calculate correlation, in pixels.
steps_cor = -size_cor:size_cor; % row vector.
Tcor = zeros(size(steps_cor,2),size(steps_cor,2));
for ii = 1:2*size_cor+1 
    dy = steps_cor(ii); % displacement along y.
    for jj = 1:2*size_cor+1 
        dx = steps_cor(jj); % displacement along x.
        cor = Bbot.*circshift(Btop,[dy,dx]); % shift the top channel by [dy dx] and multiply by bottom channel.
        Tcor(ii,jj) = sum(sum(cor)); % summ all pixels and store in result Tcor.
    end
end

% Position of maximum in the cross-correlation function gives the
% displacement (shift) of the top channel with respect to the bottom one:
[y_max x_max] = find(Tcor==max(max(Tcor))); 

% Error control:
if length(y_max)>1 || length(x_max)>1
    disp('ERROR !!!')
    disp('More than one maximum found in cross-correlation. Exiting function.')
end

% Find final shift needed in image:
% Because steps_cor is (-size_cor:size_cor), and size of cross-correlation
% function is 2*size_cor+1, we need to subtract (size_cor+1) to make the
% cross-correlation be centred at zero shift:
shift_x = x_max - (size_cor+1);
shift_y = y_max - (size_cor+1);

if show_plots ==1
    % Plot cross-correlation:
    figure;
    imshow(Tcor,[]) % show plot of cross-correlation result.
    hold on;
    plot(x_max,y_max,'x','Color','g','MarkerSize',10) % plot max of cross-correlation as a green cross on top.
    hold off;
    % Plot cross-correlation as a 3D mesh too:
    figure;
    mesh(Tcor)
end

% --------------------------------------------------
% Shift the top channel to counteract slight displacement between channels:
% the transformation that works well for images in data set 'cybD-mCherry
% &ATPase-GFp' is [5,-4]: shift down 5 pixels and left 4 pixels.
% For data set 'GFP-nuoF & mCherry-sdhC' use [3,-3] or [3,-2].
% For data set 'cydB-mCherry & GFPuv4-nuoF' use [5,-5] or [3,-3] or [4,-5].
% translateBy = [6,-1]; % amount by which to shift the image (down and left).

translateBy = [shift_y,shift_x]; % amount by which to shift the image (down and left).
Btop2 = circshift(Btop,translateBy); % shift top image down and left.

result1 = Btop2-Bbot; % shifted top minus bottom channel.
% The interesting plot is that of Btop2-Bbot, with values between -1 and 1.
% Bright areas with values close to 1 will correspond to spots found in the
% top frame but not in the bottom one. 
% Dark areas with values close to -1 will correspond to spots found in the
% bottom brame but not in the top one.
% Grey areas with values close to 0 will correspond to either the
% background or areas where bright spots are colocalised in both channels.


%% Create new directory for saving plots and results:

% Create new folder for outputs inside current folder:
output_folder_name = 'colocalisResults';
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(output_folder_name); % make new directory.


if show_plots ==1
    % Plot result (black and white):
    figure
    subplot(3,1,1); imshow(Btop2,[]); colorbar('location','eastoutside'); title('shifted top channel');
    subplot(3,1,2); imshow(Bbot,[]); colorbar('location','eastoutside'); title('bottom channel');
    subplot(3,1,3); imshow(result1,[]); colorbar('location','eastoutside'); title('shifted top minus bottom');
    
    if saveplot == 1
        % Save current figure as .png in output folder and close image:
        figName = strcat('BF_',image_label,'_frames',num2str(start_frame),'to',num2str(end_frame)); % name of figure file to save to.
        saveFigurePNG(output_folder_name,figName)
    end
end
