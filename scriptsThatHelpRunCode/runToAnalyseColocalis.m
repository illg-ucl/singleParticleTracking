% runToAnalyseColocalis.m
% cybD-mCherry-ATPase-GFp data set.
%
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%

%% First, generate plots for bright field images:
% For bright field (BF) images:

% [ybf497 xbf497] = overlapBF('497',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf499 xbf499] = overlapBF('499',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf508 xbf508] = overlapBF('508',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf512 xbf512] = overlapBF('512',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf514 xbf514] = overlapBF('514',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf517 xbf517] = overlapBF('517',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf521 xbf521] = overlapBF('521',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% [ybf523 xbf523] = overlapBF('523',1,'end',30,30,30,30,1,1);
% close all; % close all open figures
% 
% save 'shift_topChannel' 'x*' 'y*'
% 
% % Result:
% % In order of images, from the bright field images, top channel shifted by:
% [ybf497 xbf497; ybf499 xbf499; ybf508 xbf508; ybf512 xbf512; ybf514 xbf514; ybf517 xbf517; ybf521 xbf521; ybf523 xbf523]
% %      5    -4
% %      6    -5
% %      4    -7
% %      7    -2
% %      6    -5
% %      9    -7   %  translateBy = [2,-5]; could be better, obtained from  [ybf517 xbf517] = overlapBF('517',1,'end',10,20,100,100,1,1);
% %      3    -6
% %      5    -5
% -----------------------------------------


%% Colocalisation analysis of all image sequences in this data set.
% Using frame averages of both bright field image sequences and
% fluorescence image sequences.
% Reminder:
% quickGraphColocalisation(image_label,bf_image_label,start_frame,end_frame
% ,roi,output_label)

% -----------------------------------------
% For the whole image as roi ('full_image'):

% quickGraphColocalisation('498','497',7,'end','full_image','_full_image')
% quickGraphColocalisation('500','499',5,'end','full_image','_full_image')
% quickGraphColocalisation('509','508',22,'end','full_image','_full_image')
% quickGraphColocalisation('513','512',11,'end','full_image','_full_image')
% quickGraphColocalisation('515','514',10,'end','full_image','_full_image')
% quickGraphColocalisation('518','517',7,'end','full_image','_full_image')
% quickGraphColocalisation('522','521',2,'end','full_image','_full_image')
% quickGraphColocalisation('524','523',20,'end','full_image','_full_image')
% -----------------------------------------


% -----------------------------------------
% To do using a region of interest (roi), use first frameAverage to get
% coordinates of roi as [xleft xright ytop ybottom]. Use coordinates from
% the top channel (as if the image was only of total size the size of the
% top half).

% Run one by one and write the roi below:
% % frameAverage('498',7,'end',1,0); % display (1), but no saving (0).
% % frameAverage('500',5,'end',1,0);
% % frameAverage('509',22,'end',1,0);
% % frameAverage('513',11,'end',1,0);
% % frameAverage('515',10,'end',1,0);
% % frameAverage('518',7,'end',1,0);
% % frameAverage('522',2,'end',1,0);
% % frameAverage('524',20,'end',1,0);

% roi498 = [307 384 37 129];
% roi500 = [275 377 34 126];
% roi509 = [154 235 25 98];
% roi513 = [297 384 38 117];
% roi515 = [237 298 51 116];
% roi518a = [228 283 16 74];
% roi518b = [181 270 172 235];
% roi522a = [234 297 5 45];
% roi522b = [223 268 46 98];
% roi524 = [278 329 21 72];
%  

% For:
% colocalis_threshold = 0.144;
% factor_forThresholding = 1.5;

% quickGraphColocalisation('498','497',7,'end',roi498,'_roi')
% quickGraphColocalisation('500','499',5,'end',roi500,'_roi')
% quickGraphColocalisation('509','508',22,'end',roi509,'_roi')
% quickGraphColocalisation('513','512',11,'end',roi513,'_roi')  % run with translateBy = [6,-1];
% quickGraphColocalisation('515','514',10,'end',roi515,'_roi')  % run with translateBy = [3,-5];
% quickGraphColocalisation('518','517',7,'end',roi518a,'_roiA')  % run with translateBy = [2,-5];
% quickGraphColocalisation('518','517',7,'end',roi518b,'_roiB')
% quickGraphColocalisation('522','521',2,'end',roi522a,'_roiA')
% quickGraphColocalisation('522','521',2,'end',roi522b,'_roiB')
% quickGraphColocalisation('524','523',20,'end',roi524,'_roi')
% % -----------------------------------------

% cd('C:\Isabel\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\colocalisResults\ForAveragesOfAllFramesWithROI')
% mat_names_0 = dir('*.mat'); % list of all .mat file names.
% mat_names = {mat_names_0.name}'; % column array of file names.
% % Put together all colocalisation values from all frame averages of all
% % image sequences;
% all_coloc_values_allAvgImages = []; % initialise empty vector to accumulate values.
% for i=1:length(mat_names)
%     load(mat_names{i}) % this loads the variable named "colocalis_values" (same for all images) onto the workspace. (column vector)
%     all_coloc_values_allAvgImages = [all_coloc_values_allAvgImages; colocalis_values];
% end
% % plot histogram:
% hist(all_coloc_values_allAvgImages,100) % use 100 bins.
% xlabel('colocalisation value')
% ylabel('frequency')
% xlim([-0.8 0.8])
% % save:
% save 'all_colocValues_allFrameAvgsAllImages' 'all_coloc_values_allAvgImages' % save .mat file
% xlswrite('all_colocValuesFrameAvgs_allImages.xls',all_coloc_values_allAvgImages,'all_colocValuesFrameAvgs'); % save .xls file.
% % return to main path:
% cd('C:\Isabel\ExperimData\HeikoData\cybD-mCherry &ATPase-GFp\cybD-mCherry &ATPase-GFp TIRF\good\')
% 


%% Analysis frame by frame:
% eg.  rs498 = colocalisationFrameByFrame('498','497',7,7+25,7,'end',[307 384 37 129],'');
% see folder: For498first25frames_roi

% Using  colocalis_threshold = 0.144;
% Using  factor_forThresholding = 1.5;

% colocFrameByFrame498 = colocalisationFrameByFrame('498','497',7,7+24,7,'end',roi498,'_roi',0.144,1.5);
% colocFrameByFrame500 = colocalisationFrameByFrame('500','499',5,5+24,5,'end',roi500,'_roi',0.144,1.5);
% colocFrameByFrame509 = colocalisationFrameByFrame('509','508',22,22+24,22,'end',roi509,'_roi',0.144,1.5);
% colocFrameByFrame513 = colocalisationFrameByFrame('513','512',11,11+24,11,'end',roi513,'_roi',0.144,1.5);
% colocFrameByFrame515 = colocalisationFrameByFrame('515','514',10,10+24,10,'end',roi515,'_roi',0.144,1.5);
% colocFrameByFrame518A = colocalisationFrameByFrame('518','517',7,7+24,7,'end',roi518a,'_roiA',0.144,1.5);
% colocFrameByFrame518B = colocalisationFrameByFrame('518','517',7,7+24,7,'end',roi518b,'_roiB',0.144,1.5);
% colocFrameByFrame522A = colocalisationFrameByFrame('522','521',2,2+24,2,'end',roi522a,'_roiA',0.144,1.5);
% colocFrameByFrame522B = colocalisationFrameByFrame('522','521',2,2+24,2,'end',roi522b,'_roiB',0.144,1.5);
% colocFrameByFrame524 = colocalisationFrameByFrame('524','523',20,20+24,20,'end',roi524,'_roi',0.144,1.5);
%  
% save 'results_colocFrameByFrame' 'colocFrameByFrame*'

% Results in order as above:
% Left = negative = on bottom channel = green channel.
% Right = positive = on top channel = red channel.
% Centre = ~zero = on both channels.
% left_percent_fractions = [11 36 6 11 8 16 13.3 12 13 42];
% centre_percent_fractions = [6 62 0 1 0 0 41.4 7 12 6];
% right_percent_fractions = [83 2 94 88 92 84 45.3 81 75 52];
% 
% mean(left_percent_fractions) % 16.83 ~ 17
% std(left_percent_fractions)  % 12.09
% length(left_percent_fractions) % 10
% std(left_percent_fractions)/sqrt(length(left_percent_fractions)) % 3.82 ~ 4
% 
% mean(centre_percent_fractions) % 13.54 ~ 14
% std(centre_percent_fractions)  % 21.1
% length(centre_percent_fractions) % 10 
% std(centre_percent_fractions)/sqrt(length(centre_percent_fractions))  % 6.7 ~ 7
% 
% mean(right_percent_fractions)  % 69.63 ~ 70
% std(right_percent_fractions)   % 28.8
% length(right_percent_fractions) % 10
% std(right_percent_fractions)/sqrt(length(right_percent_fractions))  % 9.1 ~9

% Npoints = length(left_percent_fractions);
% x_axis = (1:Npoints);
% plot(x_axis,left_percent_fractions,'ro-')
% hold on
% plot(x_axis,centre_percent_fractions,'gs-')
% plot(x_axis,right_percent_fractions,'b^-')
% ylim([0 100])
% xlabel('image number')
% ylabel('% fraction')
% legend('left','centre','right','Location','Best')
% title('cybD-mCherry & ATPase-GFP')
% hold off
% 
% percent_fractions_matrix = [left_percent_fractions' centre_percent_fractions' right_percent_fractions'];
% % In percent_fractions_matrix, each column contains the left, centre and right percent fractions
% % respectively.
% width_histogramBars = 5;
% vector_edges_histogramBars = (0:width_histogramBars:100); % vector of edges for the histogram bars, from 0 to 100% in steps of width_histogramBars%.
% vector_centres_histogramBars = (width_histogramBars/2:width_histogramBars:100+width_histogramBars/2);
% n_freqs_matrix = histc(percent_fractions_matrix,vector_edges_histogramBars); % "n_freqs_matrix" is a matrix where each column contains the frequencies for each bin with the above edges, for left, centre and right columns.
% % function "histc" counts the number of occurrences such as it is higher or
% % equal to lower edge of bin, but lower than higher edge of bin: edges(k)
% % <= value < edges (k+1), and the very last bin counts occurrences which
% % are exactly equal to the final highest edge: value == edges(end). That's
% % why each column in n_freqs_matrix has the same length as
% % vector_edges_histogramBars.
% 
% % Plot three histograms together of the left(red), centre(green) and
% % right(blue) percentage fractions in the same figure:
% figure;
% h = bar(vector_centres_histogramBars,n_freqs_matrix);
% % Get the colours of each histogram right:
% set(h(1),'FaceColor','r')
% set(h(2),'FaceColor','g')
% set(h(3),'FaceColor','b')
% xlabel('% fraction')
% ylabel('frequency')
% set(gca,'YTick',0:5) % set vertical axis ticks to integer numbers onlys
% legend('left','centre','right','Location','Best')
% title('cybD-mCherry & ATPase-GFP')
% 

% % ------------------------
% variable = colocFrameByFrame498;
% coloc_frac = [variable.colocalised_fraction]; % row vector
% topChannel_frac = [variable.topChannel_fraction]; % row vector
% bottomChannel_frac = [variable.bottomChannel_fraction]; % row vector
% 
% length(coloc_frac)
% mean(coloc_frac)
% std(coloc_frac)
% figure; hist(coloc_frac,20); xlabel('colocalised fraction'); ylabel('frequency')
% 
% length(topChannel_frac)
% mean(topChannel_frac)
% std(topChannel_frac)
% figure; hist(topChannel_frac,20); xlabel('top channel fraction'); ylabel('frequency')
% 
% length(bottomChannel_frac)
% mean(bottomChannel_frac)
% std(bottomChannel_frac)
% figure; hist(bottomChannel_frac,20); xlabel('bottom channel fraction'); ylabel('frequency')
% % ------------------------



%% Analysis frame by frame using threshold from bgnd autoflu levels in parental strains:
% Using colocalisationFrameByFrame2.m

% Using  colocalis_threshold = 0.144; %(this doesn't matter much).
% singleFrame_topThreshold = 100;
% singleFrame_bottomThreshold = 400;

% colocFrameByFrame2_498 = colocalisationFrameByFrame2('498','497',7,7+24,7,'end',roi498,'_2_roi',0.144,1.5,100,400);
% colocFrameByFrame2_500 = colocalisationFrameByFrame2('500','499',5,5+24,5,'end',roi500,'_2_roi',0.144,1.5,100,400);
% colocFrameByFrame2_509 = colocalisationFrameByFrame2('509','508',22,22+24,22,'end',roi509,'_2_roi',0.144,1.5,100,400);
% colocFrameByFrame2_513 = colocalisationFrameByFrame2('513','512',11,11+24,11,'end',roi513,'_2_roi',0.144,1.5,100,400);
% colocFrameByFrame2_515 = colocalisationFrameByFrame2('515','514',10,10+24,10,'end',roi515,'_2_roi',0.144,1.5,100,400);
% colocFrameByFrame2_518A = colocalisationFrameByFrame2('518','517',7,7+24,7,'end',roi518a,'_2_roiA',0.144,1.5,100,400);
% colocFrameByFrame2_518B = colocalisationFrameByFrame2('518','517',7,7+24,7,'end',roi518b,'_2_roiB',0.144,1.5,100,400);
% colocFrameByFrame2_522A = colocalisationFrameByFrame2('522','521',2,2+24,2,'end',roi522a,'_2_roiA',0.144,1.5,100,400);
% colocFrameByFrame2_522B = colocalisationFrameByFrame2('522','521',2,2+24,20,'end',roi522b,'_2_roiB',0.144,1.5,100,400);
% colocFrameByFrame2_524 = colocalisationFrameByFrame2('524','523',20,20+24,20,'end',roi524,'_2_roi',0.144,1.5,100,400);
%  
% save 'results_colocFrameByFrame2_' 'colocFrameByFrame*'

% Results in order as above:
% Left = negative = on bottom channel = green channel.
% Right = positive = on top channel = red channel.
% Centre = ~zero = on both channels.
% left_percent_fractions = [];
% centre_percent_fractions = [];
% right_percent_fractions = [];
% 
% mean(left_percent_fractions) % 16.83 ~ 17
% std(left_percent_fractions)  % 12.09
% length(left_percent_fractions) % 10
% std(left_percent_fractions)/sqrt(length(left_percent_fractions)) % 3.82 ~ 4
% 
% mean(centre_percent_fractions) % 13.54 ~ 14
% std(centre_percent_fractions)  % 21.1
% length(centre_percent_fractions) % 10 
% std(centre_percent_fractions)/sqrt(length(centre_percent_fractions))  % 6.7 ~ 7
% 
% mean(right_percent_fractions)  % 69.63 ~ 70
% std(right_percent_fractions)   % 28.8
% length(right_percent_fractions) % 10
% std(right_percent_fractions)/sqrt(length(right_percent_fractions))  % 9.1 ~9

% Npoints = length(left_percent_fractions);
% x_axis = (1:Npoints);
% plot(x_axis,left_percent_fractions,'ro-')
% hold on
% plot(x_axis,centre_percent_fractions,'gs-')
% plot(x_axis,right_percent_fractions,'b^-')
% ylim([0 100])
% xlabel('image number')
% ylabel('% fraction')
% legend('left','centre','right','Location','Best')
% title('cybD-mCherry & ATPase-GFP')
% hold off
% 
% percent_fractions_matrix = [left_percent_fractions' centre_percent_fractions' right_percent_fractions'];
% % In percent_fractions_matrix, each column contains the left, centre and right percent fractions
% % respectively.
% width_histogramBars = 5;
% vector_edges_histogramBars = (0:width_histogramBars:100); % vector of edges for the histogram bars, from 0 to 100% in steps of width_histogramBars%.
% vector_centres_histogramBars = (width_histogramBars/2:width_histogramBars:100+width_histogramBars/2);
% n_freqs_matrix = histc(percent_fractions_matrix,vector_edges_histogramBars); % "n_freqs_matrix" is a matrix where each column contains the frequencies for each bin with the above edges, for left, centre and right columns.
% % function "histc" counts the number of occurrences such as it is higher or
% % equal to lower edge of bin, but lower than higher edge of bin: edges(k)
% % <= value < edges (k+1), and the very last bin counts occurrences which
% % are exactly equal to the final highest edge: value == edges(end). That's
% % why each column in n_freqs_matrix has the same length as
% % vector_edges_histogramBars.
% 
% % Plot three histograms together of the left(red), centre(green) and
% % right(blue) percentage fractions in the same figure:
% figure;
% h = bar(vector_centres_histogramBars,n_freqs_matrix);
% % Get the colours of each histogram right:
% set(h(1),'FaceColor','r')
% set(h(2),'FaceColor','g')
% set(h(3),'FaceColor','b')
% xlabel('% fraction')
% ylabel('frequency')
% set(gca,'YTick',0:5) % set vertical axis ticks to integer numbers onlys
% legend('left','centre','right','Location','Best')
% title('cybD-mCherry & ATPase-GFP')
% 

% --------------------------
% --------------------------

%% COLOCALISATION ANALYSIS TRACK BY TRACK

% analyseSetOfTracksColocalis(minNumPointsInTrack,minNumFramesOverlap)

% analyseSetOfTracksColocalis(3,3)
 
% ------------------------------
% Data set: cybD-mCherry-ATPase-GFp_
% --------
% Image sequence: 498
% There are 2 tracks with at least 3 points from top/red channel.
% There are 28 tracks with at least 3 points from bottom/green channel.
% There are 9 track-pairs which overlap in time for image sequence 498.
% --------
% Image sequence: 500
% There are 3 tracks with at least 3 points from top/red channel.
% There are 24 tracks with at least 3 points from bottom/green channel.
% There are 3 track-pairs which overlap in time for image sequence 500.
% --------
% Image sequence: 509
% There are 4 tracks with at least 3 points from top/red channel.
% There are 10 tracks with at least 3 points from bottom/green channel.
% There are 4 track-pairs which overlap in time for image sequence 509.
% --------
% Image sequence: 513
% There are 2 tracks with at least 3 points from top/red channel.
% There are 12 tracks with at least 3 points from bottom/green channel.
% There are 3 track-pairs which overlap in time for image sequence 513.
% --------
% Image sequence: 515
% There are 5 tracks with at least 3 points from top/red channel.
% There are 1 tracks with at least 3 points from bottom/green channel.
% There are 0 track-pairs which overlap in time for image sequence 515.
% --------
% Image sequence: 518
% There are 5 tracks with at least 3 points from top/red channel.
% There are 14 tracks with at least 3 points from bottom/green channel.
% There are 5 track-pairs which overlap in time for image sequence 518.
% --------
% Image sequence: 522
% There are 3 tracks with at least 3 points from top/red channel.
% There are 31 tracks with at least 3 points from bottom/green channel.
% There are 9 track-pairs which overlap in time for image sequence 522.
% --------
% Image sequence: 524
% There are 3 tracks with at least 3 points from top/red channel.
% There are 14 tracks with at least 3 points from bottom/green channel.
% There are 9 track-pairs which overlap in time for image sequence 524.
