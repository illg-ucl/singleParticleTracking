% This is a function to simulate photobleaching events in GFP

function n = bleach(N,t,tau)


% Assume that the process is markovian, Poisson distributed. 

% n_0 is the initial number of flurophores
% t is the time 
% tau is the timescale

% RandStream.setDefaultStream ...

% Create random seed:
(RandStream('mt19937ar','seed',sum(100*clock)));

n = zeros(size(t));
%nn = zeros(size(n));
%tau = 1;
% Translate into units of t./tau
T = t./tau;
R = rand(N,numel(t)); % Random number matrix, lazy method
P = exp(-T);

% For each row of R, find the element in P when P<R;

for i=1:N
    [I,J] = find(R(i,:)>=P,1,'first');
    % J is the time index when the bleaching occurs
    n(J) = n(J)+ 1;
end

n = N-[cumsum(n)];





