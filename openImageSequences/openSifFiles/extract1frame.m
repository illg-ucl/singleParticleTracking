function [frame]=extract1frame(frame_number)
%
% Created by Isabel Llorente Garcia, 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% This function chooses a .sif image sequence from an input dialog box, 
% returns the image data corresponding to its frame number "frame_number", and plots it.
% Note: only opens .sif images!
% See "extract_image_sequence_data.m" for a more general function to
% open/read image sequences.
% 
% example of how to call this function: extract1frame(100) extracts and plots
% frame 100 of the chosen image sequence and plots it.

% Data is in 'C:\Isabel\ExperimData\HeikoData\'
disp('Data is in C:\Isabel\ExperimData\HeikoData\')
disp(' ') % empty line

% uigetfile opens a file dialog box to choose image data file:
[file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
data_folder_path = strcat(path_data,file_data);
disp(data_folder_path) % write .sif image path to command window.

% first get size of .sif image file: (see IO_Input folder, SifFunctions.txt, or page 95 my notebook 1):
[ReturnCode, numFrames, ImageSize, TotalAcquisitionSize]=GetAndorSifSize(data_folder_path,0);
disp(' ') % empty line
disp(['The total number of frames in this image sequence is: ',num2str(numFrames)]) %output total number of frames to command window.
% numFrames is the length of the sequence, the number of frames.
% ImageSize is the size of the image, e.g. 512*512.
% TotalAcquisitionSize is numFrames*ImageSize

% read .sif image data. sifData is an array of 1x1 structures, as many
% columns as frames on image sequence. Reads frames 1 to numFrames:
[sifData] = read_sif_data_direct(data_folder_path,numFrames,ImageSize,1,numFrames);
% sifData is a cell with as many elements as image frames.


frame = sifData{frame_number}.sliceData; % selects the chosen frame number of the cell using {},
% and extracts frame data which is saved in the field 'sliceData'.
% % Note that the original frame up to here has intensity values which are larger than 1 and is of class single:
% min(min(frame))
% max(max(frame))
% class(frame)

% frame = mat2gray(frame);% transform matrix into a grayscale image (with values between 0 and 1) to do image operations on it.

frame = double(frame);
frame = imrotate(frame,90);
imshow(frame,[],'Border','tight');
