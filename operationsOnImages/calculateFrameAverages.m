% calculateFrameAverages.m: this is a script to calculate frame averages for the images in this folder, 
% for the corresponding frames analysed.

% frameAverage('498',7,500,1,1);
% 
% frameAverage('500',5,500,1,1);
% 
% frameAverage('509',22,500,1,1);
% 
% frameAverage('513',11,500,1,1);
% 
% frameAverage('515',10,500,1,1);
% 
% frameAverage('518',7,500,1,1);
% 
% frameAverage('522',2,500,1,1);
% 
% frameAverage('524',20,500,1,1);


% Show and save as .avi videos of moving average of image sequence
% (averaging every 3 frames):

roi498_top = [300 390 40 130]; % roi of size 90pixels by 90 pixels.
roi498_bottom = [300 390 256+40 256+130];
roi500_top = [275 377 34 126];
roi500_bottom = [275 377 256+34 256+126];
roi509_top = [154 235 25 98];
roi509_bottom = [154 235 256+25 256+98];
roi513_top = [297 384 38 117];
roi513_bottom = [297 384 256+38 256+117];
roi515_top = [237 298 51 116];
roi515_bottom = [237 298 256+51 256+116];
roi518a_top = [228 283 16 74];
roi518a_bottom = [228 283 256+16 256+74];
roi518b_top = [181 270 172 235];
roi518b_bottom = [181 270 256+172 256+235];
roi522a_top = [234 297 5 45];
roi522a_bottom = [234 297 256+5 256+45];
roi522b_top = [223 268 46 98];
roi522b_bottom = [223 268 256+46 256+98];
roi524_top = [278 329 256+21 256+72];
roi524_bottom = [278 329 256+21 256+72];

% movingAvgVideoAvi('498',7,'end',3,roi498_top,roi498_bottom,1,'') % last input =1 if you want to save .avi video.
% movingAvgVideoAvi('500',5,'end',3,roi500_top,roi500_bottom,1,'')
% movingAvgVideoAvi('509',22,'end',3,roi509_top,roi509_bottom,1,'')
% movingAvgVideoAvi('513',11,'end',3,roi513_top,roi513_bottom,1,'')
% movingAvgVideoAvi('515',10,'end',3,roi515_top,roi515_bottom,1,'')
% movingAvgVideoAvi('518',7,'end',3,roi518a_top,roi518a_bottom,1,'a')
% movingAvgVideoAvi('518',7,'end',3,roi518b_top,roi518b_bottom,1,'b')
% movingAvgVideoAvi('522',2,'end',3,roi522a_top,roi522a_bottom,1,'a')
% movingAvgVideoAvi('522',2,'end',3,roi522b_top,roi522b_bottom,1,'b')
% movingAvgVideoAvi('524',20,'end',3,roi524_top,roi524_bottom,1,'')

% movingAvgVideoAvi('498',7,'end',3,roi498_top,roi498_bottom,1,'rescaleBothChannels') % last input =1 if you want to save .avi video.
% movingAvgVideoAvi('500',5,'end',3,roi500_top,roi500_bottom,1,'rescaleBothChannels')
% movingAvgVideoAvi('509',22,'end',3,roi509_top,roi509_bottom,1,'rescaleBothChannels')
% movingAvgVideoAvi('513',11,'end',3,roi513_top,roi513_bottom,1,'rescaleBothChannels')
% movingAvgVideoAvi('515',10,'end',3,roi515_top,roi515_bottom,1,'rescaleBothChannels')
% movingAvgVideoAvi('518',7,'end',3,roi518a_top,roi518a_bottom,1,'a_rescaleBothChannels')
% movingAvgVideoAvi('518',7,'end',3,roi518b_top,roi518b_bottom,1,'b_rescaleBothChannels')
% movingAvgVideoAvi('522',2,'end',3,roi522a_top,roi522a_bottom,1,'a_rescaleBothChannels')
% movingAvgVideoAvi('522',2,'end',3,roi522b_top,roi522b_bottom,1,'b_rescaleBothChannels')
% movingAvgVideoAvi('524',20,'end',3,roi524_top,roi524_bottom,1,'rescaleBothChannels')

% roi498_top2 = [325 367 52 108];
% roi498_bottom2 = [323 365 312 368];
% movingAvgVideoAvi('498',7,'end',3,roi498_top2,roi498_bottom2,1,'_2')
% movingAvgVideoAvi('498',125,375,3,roi498_top2,roi498_bottom2,1,'_3_5to15sec')
% movingAvgVideoAvi('498',10,300,3,roi498_top2,roi498_bottom2,1,'_frames10to300')


