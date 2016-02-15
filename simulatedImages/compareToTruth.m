function results = compareToTruth(image_label,real_Xcentre,real_Ycentre,real_SigmaSpot,real_offsetBgNoise,real_bgNoiseStd,real_num_molecules)
%
% Created by Isabel Llorente-Garcia, July 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Compare the results from bright spot detection by the software on
% simulated images with the real parameters (spot position, width,
% brightness, noise, etc.) entered in function simulateImageSequence3.m.
%
% Open the excel file which contains all found trajectories (file ending in
% fullTrajs.xls) for the corresponding image sequence (image_label) and get
% the mean % error differences with the true values from the simulated image
% sequence, as set in function simulateImageSequence3.m.
% 
% INPUTS:
% - image_label: labels image sequence and excel file of results. Eg: '11'.
% - real_Xcentre, real_Ycentre: real X and Y position of the single spot in
% the simulated image, in pixels. Eg: 25, 25.
% - real_SigmaSpot: real width (sigma) of the bright spot in the simulated
% image sequence, in pixels. Eg: 3.5.
% - real_offsetBgNoise: offset background noise (0 generally for simulated
% images).
% - real_bgNoiseStd: real background noise standard deviation.
% - real_num_molecules: real number of molecules per bright spot. The real
% brightness (I0_fit, amplitude of Gaussian) of the bright spot in the
% simulated image sequence is equal to the number of molecules per spot, since, in function simulateImageSequence3.m,
% the brightness of a single molecule is I0 = 1, and when a single bright
% spot contains more than one molecule, the brightness of the spot
% increases linearly. Therefore I0_fit is essentially the number of molecules per bright spot.
%
% OUTPUTS:
% results is a structure with fields
%     'localis_errorX_pix_abs'
%     'localis_errorY_pix_abs'
%     'localis_errorX_percent_abs'
%     'localis_errorY_percent_abs'
%     'Sigma_error_percent_abs'
%     'IspTot_error_percent_abs'
%     'NumMolecs_error_percent_abs'
%     'BgNoiseStd_error_percent_abs'
%     'SNR_error_percent_abs'
%     'localis_errorX_pix'
%     'localis_errorY_pix'
%     'localis_errorX_percent'
%     'localis_errorY_percent'
%     'Sigma_error_percent'
%     'IspTot_error_percent'
%     'NumMolecs_error_percent'
%     'BgNoiseStd_error_percent'
%     'SNR_error_percent'
% 
% The first block is for mean errors (mean of absolute value of error).
% The second block is without taking absolute value, to detect systematic
% over or underestimations...
%
% Example of how to call this function:
% r03 = compareToTruth('03',25,27,3.5,0.2,1).


% Find name of the analysis excel file corresponding to the simulated image
% sequence given by "image_label":
xlsFileName0 = dir(strcat('*',image_label,'fullTrajs.xls')); % Get names of excel files in that folder (accepted track-analysis .xls files).
xlsFileName = xlsFileName0.name;

% The cell array "sheet_names" contains the names of the sheets within the excel
% file (this is needed because if there is only one track found in the
% whole image, the results get sent to a Sheet named 'Sheet1' instead of 'Track results').
[typ, sheet_names] = xlsfinfo(xlsFileName);
% Now check if there is a Sheet named 'Track results' in the excel file
% (ismember('Track results', sheet_names) will be 1 if yes and 0 if not):
if ismember('Track results', sheet_names) == 1
    correct_sheet_name = 'Track results';
else
    correct_sheet_name = 'Sheet1';
end


%% Import the data (column vectors) in the sheet named 'Track data':
[numeric1,txt1,raw1] = xlsread(xlsFileName,correct_sheet_name);
% Turn imported data from excel file into a structure:
str_TrackData = struct; % create empty structure to fill up.
for u = 1:length(txt1)
    str_TrackData = setfield(str_TrackData, txt1{u}, numeric1(:,u)); % create field in the structure.
    % Each field in the structure contains a column vector with the data.
end
% Import all values as column vectors:
% FrameNumber = str_TrackData.FrameNumber; % column vector.
CentreX = str_TrackData.CentreX; % column vector.
CentreY = str_TrackData.CentreY; % column vector.
IspTot = str_TrackData.IspTot; % column vector.
I0Fit = str_TrackData.I0Fit; % column vector.
SigmaFit = str_TrackData.SigmaFit; % column vector.
OffsetBgNoise = str_TrackData.bg_noise_offset_afterBGsubtract;
BgNoiseStd = str_TrackData.BgNoiseStd; % column vector.
SNR = str_TrackData.SNR; % column vector.
                                
% Real values of parameters in simulated image are:
% real_Xcentre, real_Ycentre, real_SigmaSpot, real_bgNoiseStd,
% real_num_molecules.


%% CALCULATE MEAN ERRORS:

% Localisation error in pixels for X and Y:
results.localis_errorX_pix_abs = mean(abs(CentreX-real_Xcentre)); 
results.localis_errorY_pix_abs = mean(abs(CentreY-real_Ycentre)); 
% Localisation error as a percentage error:
results.localis_errorX_percent_abs = mean(100*abs(CentreX-real_Xcentre)/real_Xcentre);
results.localis_errorY_percent_abs = mean(100*abs(CentreY-real_Ycentre)/real_Ycentre);

% Spot width error as percentage error:
results.Sigma_error_percent_abs = mean(100*abs(SigmaFit-real_SigmaSpot)/real_SigmaSpot);

% Integrated spot intensity error (area under 2D Gaussian is
% 2*pi*I0*sigma^2):
real_IspTot = 2*pi*real_num_molecules*(real_SigmaSpot)^2; % noiseless.
results.IspTot_error_percent_abs = mean(100*abs(IspTot-real_IspTot)/real_IspTot);

% Number of molecules per bright spot percentage error:
% Note: the Gaussian intensity amplitude of a single molecule is 1. And if
% there are more than one molecule per bright spot, the total instensity
% amplitude (I0) of the spot is then equal to real_num_molecules.
results.NumMolecs_error_percent_abs = mean(100*abs(I0Fit-real_num_molecules)/real_num_molecules);

% Offset background noise error, absolute value, absolute error:
results.OffsetBgNoise_error_abs = mean(abs(OffsetBgNoise-real_offsetBgNoise)); % cannot do percent error comparing to zero: cannot divide by 0.

% Background noise standard deviation error as percentage error:
results.BgNoiseStd_error_percent_abs = mean(100*abs(BgNoiseStd-real_bgNoiseStd)/real_bgNoiseStd);

% Signal-to-Noise-Ratio error as percentage error:
real_SNR = real_num_molecules/real_bgNoiseStd; % I0 ampli of whole bright spot is equal to num_molecules (I0 ampli of a single molecule is 1).
results.SNR_error_percent_abs = mean(100*abs(SNR-real_SNR)/real_SNR);



% And same percentage errors but with no absolute value in their
% calculation, to detect if we systematically over or underestimate a given
% parameter:
results.localis_errorX_pix = mean(CentreX-real_Xcentre); 
results.localis_errorY_pix = mean(CentreY-real_Ycentre); 

results.localis_errorX_percent = mean(100*(CentreX-real_Xcentre)/real_Xcentre);
results.localis_errorY_percent = mean(100*(CentreY-real_Ycentre)/real_Ycentre);

results.Sigma_error_percent = mean(100*(SigmaFit-real_SigmaSpot)/real_SigmaSpot);

results.IspTot_error_percent = mean(100*(IspTot-real_IspTot)/real_IspTot);

results.NumMolecs_error_percent = mean(100*(I0Fit-real_num_molecules)/real_num_molecules);

results.OffsetBgNoise_error = mean(OffsetBgNoise-real_offsetBgNoise); % cannot do percent error comparing to zero: cannot divide by 0.

results.BgNoiseStd_error_percent = mean(100*(BgNoiseStd-real_bgNoiseStd)/real_bgNoiseStd);

results.SNR_error_percent = mean(100*(SNR-real_SNR)/real_SNR);