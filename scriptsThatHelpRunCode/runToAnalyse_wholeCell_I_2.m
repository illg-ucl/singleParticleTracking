% runToAnalyse_wholeCell_I.m
%
% Isabel, Junio 2012.
%
% Find exponential photobleaching time of the intensity integrated over the
% whole cell:

% Using this function:
% results = fullCellIntensity(dataSet_label,image_label,start_frame,tsamp,roi,output_label) 
% --------------------------------


% For GFP, GREEN channel, bottom half of image:

w498_green = fullCellIntensity('cybD-mCherry-ATPase-GFp','498',7,0.04,[1 512 257 512],'_green');

% For image 500 there are two cells too close together.

w509_green = fullCellIntensity('cybD-mCherry-ATPase-GFp','509',22,0.04,[1 512 257 512],'_green');

% For image 513 there are two cells too close together.

w515_green = fullCellIntensity('cybD-mCherry-ATPase-GFp','515',10,0.04,[1 512 257 512],'_green');

w518_green = fullCellIntensity('cybD-mCherry-ATPase-GFp','518',7,0.04,[1 512 257 512],'_green');

w522_green_cell1 = fullCellIntensity('cybD-mCherry-ATPase-GFp','522',2,0.04,[230 288 264 303],'_green_cell1');
w522_green_cell2 = fullCellIntensity('cybD-mCherry-ATPase-GFp','522',2,0.04,[220 265 308 360],'_green_cell2');

w524_green = fullCellIntensity('cybD-mCherry-ATPase-GFp','524',20,0.04,[1 512 257 512],'_green');
% --------------------------------

% For mCherry, RED channel, bottom half of image:

w498_red = fullCellIntensity('cybD-mCherry-ATPase-GFp','498',7,0.04,[1 512 1 250],'_red');

% For image 500 there are two cells too close together.

w509_red = fullCellIntensity('cybD-mCherry-ATPase-GFp','509',22,0.04,[1 512 1 250],'_red');

% For image 513 there are two cells too close together.

w515_red = fullCellIntensity('cybD-mCherry-ATPase-GFp','515',10,0.04,[1 512 1 250],'_red');

w518_red = fullCellIntensity('cybD-mCherry-ATPase-GFp','518',7,0.04,[1 512 1 120],'_red');

w522_red_cell1 = fullCellIntensity('cybD-mCherry-ATPase-GFp','522',2,0.04,[237 291 6 43],'_red_cell1');
w522_red_cell2 = fullCellIntensity('cybD-mCherry-ATPase-GFp','522',2,0.04,[227 270 49 95],'_red_cell2');

w524_red = fullCellIntensity('cybD-mCherry-ATPase-GFp','524',20,0.04,[1 512 1 250],'_red');
% --------------------------------


save 'results_wholeCell_I' 'w*' % save all result structures in a .mat file.