function yVectorsToPlot = getErrorVectors(results_vector,numMolecs_vector)
%
% Created by Isabel Llorente-Garcia, July 2012.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% Group results and plot mean errors.
%
% INPUTS: 
% - results_vector: vector such as [r11 r12 r13], where r11 is the output
% structure for a given test, output of function "compareToTruth.m" (eg. r11 = compareToTruth('11',25,25,3.5,0.5,1);).
% - numMolecs_vector: column vector with number of molecules as it was varied for
% results_vector. Eg. [1 2 5]' (1 molec for r11, 2 for r12, 3 for r13).
% 
    
field_names = fields(results_vector(1)); % cell structure

for i=1:length(results_vector)
    
    for j=1:length(field_names)
    
        yVectorsToPlot(i,j) = getfield(results_vector(i),field_names{j}); % each column corresponds to a variable j.

    end
end

% plot(numMolecs_vector,yVectorsToPlot(:,1))

