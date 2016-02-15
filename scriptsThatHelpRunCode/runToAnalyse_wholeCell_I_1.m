% runToAnalyse_wholeCell_I.m
%
% Isabel, Nov 2012.
%
% Find exponential photobleaching time of the intensity integrated over the
% whole cell:

% Using this function:
% results = fullCellIntensity(dataSet_label,image_label,start_frame,tsamp,roi,output_label) 
% --------------------------------

% Use frame averages as help to determine region of interest (roi), e.g.:
% A = frameAverage('3ob',1,'end',1,0);

% --------------------------------

% For mCherry, RED channel, whole image (this is new Oxphos data taken on
% October 24th 2012, image is not split in top and bottom channels, the
% full image is red channel):

% w3_red = fullCellIntensity('cyoA-mCherry','3ob',11,0.04,[168 232 143 208],'_red');
% 
% w5_red = fullCellIntensity('cyoA-mCherry','5ob',11,0.04,[176 221 170 247],'_red');
% 
% w6_red = fullCellIntensity('cyoA-mCherry','6ob',8,0.04,[181 261 235 305],'_red');
% 
% w8_red = fullCellIntensity('cyoA-mCherry','8ob',11,0.04,[74 150 190 256],'_red');
% 
% w10_red = fullCellIntensity('cyoA-mCherry','10ob',14,0.04,[196 243 146 206],'_red');
% % use only the non-wobbly cell in image '10ob'.
% 
% save 'results_wholeCell_I' 'w*' % save all result structures in a .mat file.

% --------------------------------
%%  Fit results for all images:

% For exponential fits of intensity vs time with no offset:
tau_fits_noOffset = [1.58 2.42 3.03 1.76 1.84]; % tau in seconds.
Error_tau_fits_noOffset = [0.03 0.04 0.03 0.03 0.03]; % error of stau in seconds.

mean(tau_fits_noOffset) % 2.1260
std(tau_fits_noOffset) % 0.5953
0.5953/sqrt(5) % 0.2662

sqrt(sum(Error_tau_fits_noOffset.^2)) % 0.0721


% For exponential fits of intensity vs time with offset:
tau_fits_offset = [1.20 1.63 2.25 1.23 1.33]; % tau in seconds.
Error_tau_fits_offset = [0.01 0.02 0.02 0.013 0.02]; % error of stau in seconds.

mean(tau_fits_offset) % 1.5280
std(tau_fits_offset) % 0.4380
0.4380/sqrt(5) % 0.1959

sqrt(sum(Error_tau_fits_offset.^2)) % 0.0383

% -----------------


% Use function “analyseSetOfExcelFiles()” to produce the table of results:

% analyseSetOfExcelFiles()

% Now back fit all track intensities using a fixed tau of 2.1 plus minus 0.6 seconds:
% Use a single molecule in vitro mCherry intensity of 500 +- 600 (the error bar (600) is not used for calculations).

% analyseSetOfTracks2(minNumPointsInTrack,max_framesAwayFromTimeOrigin,colo
% ur_chosen,tau_fixed,min_rsq_fit_I,max_fit_error,prob_accept_fit,I_singleMolec,error_IsingleMolec)

% No limits first;
analyseSetOfTracks2(3,500,'any',2.1,-Inf,Inf,0,500,600) % 20 tracks.
% The number of tracks fulfilling conditions is: 28.
% See folder "cyoA_mCherry_I0_any_noLimits"

% Now introducing constraints:
analyseSetOfTracks2(3,500,'any',2.1,0,Inf,0.01,500,600)
% The number of tracks fulfilling conditions is: 14.
% See folder "cyoA_mCherry_I0_any_480pm50"


