%% DEFINITIONS and PARAMETERS:

% file "paramsForLinkTrajSegments.m"

% Use a very large number (larger than image size in pixels) for rejected asignments:
rej = 10000;
params.rej = rej; % Save parameters to results as structure "params".
 
% PARAMETERS for linking trajectory segments:
% For linking end spot in one trajectory with start spot in another trajectory:
d_01_max = 5; % max distance in pixels between spot centres.(5)
Iratio_01_min = 0.2; % min ratio of total spot intensities (after bgnd subtraction).(0.2)
Iratio_01_max = 4; % max ratio of total spot intensities (after bgnd subtraction). (Large enough to account for blinking, note that Iratio is for frame k-1/frame k.)(3)
SigmaRatio_01_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit).(0.5)
SigmaRatio_01_max = 2; % max ratio of spot width (sigma of Gaussian fit). (2)
Frames_away_max = 2; % max separation in frames (i.e., prop to time) for trajectory segments to be linked.(default = 2). Write 1 if you want no jumps (no frame jumps) in the segment linking.
% At most we skip one frame when Frames_away_max = 2.
% Save parameters to results as structure "params":
params.d_01_max = d_01_max; 
params.Iratio_01_min = Iratio_01_min; 
params.Iratio_01_max = Iratio_01_max; 
params.SigmaRatio_01_min = SigmaRatio_01_min;
params.SigmaRatio_01_max = SigmaRatio_01_max; 
params.Frames_away_max = Frames_away_max;
