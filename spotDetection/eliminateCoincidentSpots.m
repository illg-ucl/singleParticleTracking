function [new_spotsX new_spotsY candid_pos_to_keep] = eliminateCoincidentSpots(spotsX,spotsY,limit_dist)
%
% Created by Isabel Llorente-Garcia, 2011.
% If you use this code please acknowledge Isabel Llorente-Garcia in your
% publications.
%
% ELIMINATE COINCIDENT points
%
% Inputs: spotsX and spotsY are two column vectors with the x and y
% coordinates of the spots, respectively.
% 
% The function checks the distances between all pairs of points and removes
% those points (x,y) which are closer than limit_dist pixels (distance<limit_dist) to another point in the list. 


% Put together in a matrix with two columns with x and y positions of
% all candidates:
allCandidates = [spotsX spotsY];
P = size(allCandidates,1);

% Error control: need at least two spots to compare:
if P<2
    new_spotsX = spotsX;
    new_spotsY = spotsY;
    if P==1
        candid_pos_to_keep = 1;
    elseif P==0
        candid_pos_to_keep = [];
    end
    return
end

distCandidates = pdist(allCandidates); % Calculate Euclidean distance between pairs of candidate spots.
link_result = linkage(distCandidates); % Link candidate spots according to their distances.
% link_result is a matrix where the first two columns identify the spots
% that have been linked and the third column contains distance between them.

% Eliminate coincident spot candidates which are closer than one pixel to each other. 
% Find pairs of candidates with distances above or equal to 1 pixel and keep them all (see to_keep1). 
% Find pairs of candidates with distances below 1 pixel and keep only one of the pair (see to_keep2).
% Find solitary spots with distances above or equal to 1 pixel to any other spot (see to_keep3).
to_keep1 = find(link_result(:,3)>=limit_dist & link_result(:,1)<=P & link_result(:,2)<=P);
to_keep2 = find(link_result(:,3)<limit_dist & link_result(:,1)<=P & link_result(:,2)<=P);
to_keep3 = find(link_result(:,3)>=limit_dist & xor(link_result(:,1)<=P,link_result(:,2)<=P));
intermediate = link_result(to_keep3,[1 2]);

% For debugging:
% link_result(to_keep1,1)
% link_result(to_keep1,2)
% link_result(to_keep2,1)
% intermediate(intermediate<=P)

% candidates positions to be kept:
candid_pos_to_keep = sort([link_result(to_keep1,1); link_result(to_keep1,2); link_result(to_keep2,1); intermediate(intermediate<=P)]);

% Candidates we keep:
new_spotsX = spotsX(candid_pos_to_keep);
new_spotsY = spotsY(candid_pos_to_keep);

% disp(['no. of initial spots in list: ',num2str(length(spotsX))])
% disp(['no. of final spots after eliminating coincidences in list: ',num2str(length(new_spotsX))])