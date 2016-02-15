function [I0vector I0stdVector] = extractI0(processedManyTraj)
%
% Created by Isabel Llorente-Garcia, Nov 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% 
% This function selects the results which have a 'GoodTrajFlag' equal to 1,
% and puts together all the I0 (intensity at time origin) values in a
% vector 'IOvector', which is the output. 
%
% Input: 'processedManyTraj' is a cell array. It is the output of function
% "showManyTrajAnalysis.m" and contains all the trajectory analysis results 
% for a number of trajectories in a given image sequence.

% Initialise output as empty vector:
I0vector = [];
I0stdVector = [];

% Loop through results for each trajectory i:
for i = 1:length(processedManyTraj) 
    flag = processedManyTraj{i}{1}.GoodTrajFlag;
    
    if flag == 1
        I0_value = processedManyTraj{i}{1}.IntensityAtTimeOrigin;
        I0_std = processedManyTraj{i}{1}.StDevI0;
        I0vector = [I0vector I0_value];
        I0stdVector = [I0stdVector I0_std];
    end
end

% Save result (as .mat) in a folder specified by user:
uisave({'I0vector','I0stdVector'},'I0results')