
function tracks = getDisplacement(tracks)
%
% Function to find how the particle moves - i.e. it is the deviation (not
% squared)
%
% Created by Alex Robson.
% Commented and modified by Isabel 18 October 2011.
%
% The input tracks is the output of function parseXLS2mat.
% The output returns the trajectories including the msd, delta t, etc...

% -----------------------------------------------
% PARAMETERS: 

tsamp = 1; % Set this to the correct sampling time. This is used to calibrate the absolute time later on, to go from frames to time in seconds. It is the time between frames.

% -----------------------------------------------

N = tracks.numel; % Number of data points in a track.
Ndiff = nchoosek(N,2); % Ndiff is the number of combinations of N things taken 2 at a time.

% Get the pairwise difference matrix:
x = tracks.xvalues;
y = tracks.yvalues;
t = double(tracks.timeabs); % Note (Isabel): use the absolute time.

PDMx = getPDM(x); % GET Pairwise Difference Matrix. See below, function getPDM is at the end of this file. 
PDMy = getPDM(y);
PDMt = getPDM(t);

% Place the non-zero values into a vector:
tdiff = PDMt(find(PDMt~=0));
ydiff = PDMy(find(PDMy~=0));
xdiff = PDMx(find(PDMx~=0));

[I,J] = sort(tdiff,'ascend');
% I is the value, J is the index. 
% The displacements over each time interval actually correspond to each diagonal in getPDM.

ydiff = ydiff(J(:));
tdiff = single(tdiff(J(:)));
xdiff = xdiff(J(:));

% now find unique values of tdiff (Delta t values).
[L,M] = unique(single(tdiff),'first');
% [L,M] = UNIQUE(X,'first') returns the unique values of vector X in vector L, and output vector M gives indices for the first occurrence of each unique value in X.  
% L has no repetitions and is also sorted.

% Check that there are no more than N-1 elements in M. Required because
% numerical accuracy means unique was failing:
if numel(M)>Ndiff 
    error('getDisplacement.m failed. Check the numerical accuracy of the calculation'); 
end

% Now place each unique value of t into a displacement vector. 
% Get number of unique values
numunq = numel(M); % number of different Delta t values.

for k=1:(numunq-1) % exclude the last point which looks only at first and last point from the analysis.

    disp{k}.xdiff = xdiff(M(k):M(k+1)-1);
    disp{k}.ydiff = ydiff(M(k):M(k+1)-1);
    disp{k}.tdiff = tdiff(M(k));
   
    % Calculate the average MSD:
    msd(k) = mean(disp{k}.xdiff.^2+disp{k}.ydiff.^2);
    deltaTime(k) = disp{k}.tdiff;
    % Relative error of msd (see "Single particle tracking". H. Quian et
    % al. Biophys. J. 60, 910-921, 1991). N is the number of data points in
    % the trajectory and n = deltaTime(k)/tsamp is the deltaTime in frames:
    n = deltaTime(k)/tsamp;
    errorRelMsd(k) = sqrt((2*n^2+1)/(3*n*(N-n+1)));
    % Note that the previous relative error assumes that all points in the
    % original trajectory are equally spaced in time/frames... So it is
    % only an approximation if we have unregularly spaced points in a
    % trajectory.
    % Absolute error of msd:
    errorMsd(k) = msd(k)*errorRelMsd(k);
end

tracks.timediff = deltaTime';
tracks.msd = msd';
tracks.errorMsd = errorMsd';
tracks.errorMsdRelPercent = 100.*errorRelMsd';
tracks.disp = disp;

% This is the original code - assuming uniformly sampled tracks (constant frame spacing). 
% % Get MSD2
% 
% for i=1:length(tracks)
%   temp2 = zeros(length(tracks(i).time),1);
%     for n=0:length(tracks(i).time)-1 % timei ndexes
%       
%     N = length(tracks(i).time);
% 
%     for j=1:N-n
%        % tracks.disp2{n+1}(j,:) = [tracks(i).xvalues(j+n) tracks(i).xvalues(j) tracks(i).xvalues(j+n)-tracks(i).xvalues(j)];
%        temp2(n+1)=temp2(n+1)+(tracks(i).xvalues(j+n)-tracks(i).xvalues(j)).^2+(tracks(i).yvalues(j+n)-tracks(i).yvalues(j)).^2;
%     end
%     temp2(n+1) = temp2(n+1)/(N-n);
%     end
% 
%     msd = temp2;     
% end 
% 
% tracks.msd2= msd;
% 
% end

end
    
    
function PDM = getPDM(x)
    % Get pairwise difference matrix:
PDM = tril(repmat(x,[1 length(x)])-repmat(x',[length(x),1]),-1); 
% TRIL(X,K) is the elements on and below the K-th diagonal
% of X .  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
% REPMAT(A,[M N]) creates a large matrix consisting of an M-by-N
% tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
end

% Examples for getPDM:
% x = [1 3 5 6]'; % column vector.
% repmat(x,[1 length(x)]) =
%      1     1     1     1
%      3     3     3     3
%      5     5     5     5
%      6     6     6     6
% repmat(x',[length(x),1]) =
%      1     3     5     6
%      1     3     5     6
%      1     3     5     6
%      1     3     5     6
% repmat(x,[1 length(x)])-repmat(x',[length(x),1]) =
%      0    -2    -4    -5
%      2     0    -2    -3
%      4     2     0    -1
%      5     3     1     0
% getPDM(x) = 
%      0     0     0     0
%      2     0     0     0
%      4     2     0     0
%      5     3     1     0
% The displacements over each time interval actually correspond to each
% diagonal in getPDM.






