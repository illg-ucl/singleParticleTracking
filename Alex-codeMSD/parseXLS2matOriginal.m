% This is a file to parse the structures created by Isabels tracking code
% (29/09/2011) to a struct file as used by Alex

% Revised - 14-10-11

function tracks = parseXLS2mat(file)


% Assumptions
tsamp = 1; % Set this to the correct sampling time. 

% Dependencies: xlsread.mat
% Check dependendies:
if 2~=exist('xlsread') error('Check file dependencies - you need to install xlsread'); end

% Open the data

[NUMERIC,TXT,RAW]=xlsread(file);

% Check that the column tags are organized as in this version

% Import the column heads and assign ID
colheads = TXT;

% Assign ID to X,Y,Int and Time.

%colheads = {'CentreX','CentreY','IspTot','rsqFit','SigmaFit','I0Fit','BgNoiseStd','IbgAvg','IbgTot','SNR','IinnerTot','ClipFlag','TrajNumber','FrameNumber','SpotNumber'};
% Generate ID
for i=1:numel(colheads) 
    ID.(colheads{i}) = find(strcmp(TXT,colheads{i})); 
end






% All comparisons should have a strcmp of 1, i.e. all not equal should be
% false, so any should be false.

if true(any(~strcmp(colheads,TXT)))     error('Incorrect mapping between column heads read and expected inputs. Check the column arrangement!'); end

% Write into a structure file as used by Alex

% Do some error checking on the file

% The trajectory column SHOULD produce incremental trajectory outputs.
% Check this.

traj = NUMERIC(:,ID.TrajNumber);
temp = diff(NUMERIC(:,ID.TrajNumber));
% Remove zero elements
temp(temp==0)=[];

% This should be an logically true array. Check this.
if ~any(logical(temp)) error('Input parsing failed. Please check the output of the trajectory file. It is a requirement that the trajectories are incrementally ordered'); end

disp('File Loaded successfully!');

% Get individual tracks
% List the points at which the trajectory first appears
[A,I,J] = UNIQUE(traj,'first');
[B,Y,Z] = UNIQUE(traj,'last');
% Get the number of tracks
numtracks = numel(B);

% Delete tracks that are less than 5 data points long.


% Prepare structure
%tracks = struct(numtracks);
a = 1;

% Create tracks structure

tracks(1:numtracks) = struct('XLS',[],'xvalues',[],'yvalues',[],'intensity',[],'msd_unavg',[],'timeabs',[],'time',[],'numel',[],'disp',[],'msd',[],'timediff',[]);

for i=1:numtracks 
    a = I(i);
    b = Y(i); 
    
    if b-a>=5 % only analyse tracks which have at least 5 points (frames) in them.
    
    data{i} = NUMERIC(a:b,:);
    tracks(i).XLS.track_index = B(i)
    tracks(i).xvalues = data{i}(2:end,ID.CentreX);
    tracks(i).yvalues = data{i}(2:end,ID.CentreY);
    % Set origin to zero
    tracks(i).xvalues = tracks(i).xvalues - (tracks(i).xvalues(1));
    tracks(i).yvalues = tracks(i).yvalues - (tracks(i).yvalues(1));
    tracks(i).intensity = data{i}(2:end,ID.I0Fit);
    tracks(i).msd_unavg = tracks(i).xvalues.^2+tracks(i).yvalues.^2;
    tracks(i).timeabs = data{i}(2:end,ID.FrameNumber).*tsamp;
    tracks(i).time = tracks(i).timeabs-tracks(i).timeabs(1);
    tracks(i).numel = b-a;
    tracks(i) = getDisplacement(tracks(i));
    else
%         tracks(i) = [];
        % delete index
        DEL(i) = i;
    end
    
end

tracks(find(DEL))=[];



% Examples of how to run this function:
% Go to the directory where the .xls file which contains the trajectory data is and then for example, make a FILELIST:
% FILELIST = dir('*.xls')
% Then to analyse the first excel file in that list do:
% output = parseXLS2mat(FILELIST(1).name);
% output is a structure with fields XLS
%     XLS: (can obtain traj number from here)
%     xvalues 
%     yvalues
%     intensity
%     msd_unavg 
%     timeabs: absolute time (equal to frame number here)
%     time (ignore this)
%     numel : number of points in the track
%     disp : complete list of pair-wise displacements before averaging for
%     each delta t.
%     msd: average mean square displacements for each delta t
%     timediff: delta t: time difference (in frames)
% Example of how to produce plots:    
% subplot(2,2,1);plot(output(2).timeabs,output(2).xvalues);title('x vs frame number');...
% subplot(2,2,2);plot(output(2).timeabs,output(2).yvalues);title('y vs frame number');...
% subplot(2,2,3);plot(output(2).timeabs,output(2).intensity,'g--x');title('Intensity vs frame number');...
% subplot(2,2,4);plot(output(2).timediff,output(2).msd,'r--x');title('Avge Mean Square Displacements vs Delta t (frames)');
% 
% 
