
function tracks = getDisplacement(tracks,tsamp)
%
% Calculate mean square displacements as a function of time-interval delta t.
%
% Created by Alex Robson.
% Commented and modified by Isabel October 2011.
% If you use this code please acknowledge Alex Robson and Isabel Llorente-Garcia in your
% publications.
%
%
% The input 'tracks' is an intermediate output of function 'analyseTraj'.
% 'tsamp' is the sampling time, used to calibrate the absolute time, to go from frames to time in seconds. 
% It is the time between frames in seconds. 
% The output returns the trajectory data including the msd, delta t, msd
% error, etc...

N = tracks.numel; % Number of data points in a track.
Ndiff = nchoosek(N,2); % Ndiff is the number of combinations of N things taken 2 at a time.

% Get the pairwise difference matrix:
x = tracks.xvalues;
y = tracks.yvalues;
t = tracks.timeabs; % Note (Isabel): use the absolute time.

PDMx = getPDM(x); % GET Pairwise Difference Matrix. See getPDM.m. 
PDMy = getPDM(y);
PDMt = getPDM(t);

t1 =ones(N,N);
t2 = logical(tril(t1,-1)); % triangular matrix with logical ones in lower -1 triangle. Marks positions of numbers to take from the previous matrices.

% Place the non-zero values into a vector:
tdiff = PDMt(t2);
ydiff = PDMy(t2);
xdiff = PDMx(t2);
% % tdiff = PDMt(find(PDMt~=0));
% % ydiff = PDMy(find(PDMy~=0));
% % xdiff = PDMx(find(PDMx~=0));

[I,J] = sort(tdiff,'ascend');
% I is the value, J is the index. 
% The displacements over each time interval actually correspond to each diagonal in getPDM.

ydiff = ydiff(J(:));
tdiff = single(tdiff(J(:)));
xdiff = xdiff(J(:));

% now find unique values of tdiff (Delta t values).
[L,M] = unique(single(tdiff),'first');
% [L,M] = UNIQUE(X,'first') returns the unique values of vector X in vector L, 
% and output vector M gives indices for the first occurrence of each unique value in X.  
% L has no repetitions and is also sorted.

% Check that there are no more than N-1 elements in M. Required because
% numerical accuracy means unique was failing:
if numel(M)>Ndiff 
    error('getDisplacement.m failed. Check the numerical accuracy of the calculation'); 
end

% Now place each unique value of t into a displacement vector. 
% Get number of unique values
numunq = numel(M); % number of different Delta-t values.

if numunq >1
    for k=1:(numunq-1) % exclude the last point which looks only at first and last points from the analysis.
        
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
        % Absolute error of msd (sqrt of variance, stdev):
        errorMsd(k) = msd(k)*errorRelMsd(k);
    end
else % if there is only one value of deltaTime (ie if track only has two points):
    disp{1}.xdiff = xdiff;
    disp{1}.ydiff = ydiff;
    disp{1}.tdiff = tdiff(1);
    
    % Calculate the average MSD:
    msd(1) = mean(disp{1}.xdiff.^2+disp{1}.ydiff.^2);
    deltaTime(1) = disp{1}.tdiff;
    % Relative error of msd (see "Single particle tracking". H. Quian et
    % al. Biophys. J. 60, 910-921, 1991). N is the number of data points in
    % the trajectory and n = deltaTime(k)/tsamp is the deltaTime in frames:
    n = deltaTime(1)/tsamp;
    errorRelMsd(1) = sqrt((2*n^2+1)/(3*n*(N-n+1)));
    % Note that the previous relative error assumes that all points in the
    % original trajectory are equally spaced in time/frames... So it is
    % only an approximation if we have unregularly spaced points in a
    % trajectory.
    % Absolute error of msd (sqrt of variance, stdev):
    errorMsd(1) = msd(1)*errorRelMsd(1);
end

tracks.deltaTime = deltaTime';
tracks.msd = single(msd'); % convert to single precision.
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
    
 

