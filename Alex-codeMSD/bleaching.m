% Function that models bleaching effects
function bleaching


T = 20; % Number of tracks that are being sampled - make this a square number!
N = 10; % Number of protomers in one protein complex
T_max = 1;
T_res = 0.004;
tau = 0.8;
t_rat = linspace(0,T_max/tau,round(T_max/T_res));
S = round(T_max/T_res); % Number of samples in each track
I_1 = 100;

add_noise=1; % This parameters adds in noise to the intensity of each protomer
add_bg_noise=1; % This parameter adds in background noise (requires add_noise=1)
bg = 1000; % Intensity of the background noise
pn = 100; % Scale factor of the protomor noise



for t=1:T
        r = rand(N,1);
    for i=1:N
        P = exp(-t_rat); %% MAKE THIS POISSON
        
        time = find(P<r(i),1,'first');
        if ~isnan(time)
        bleaching_time_index(i) = time; % This is the index of when the protomer is photobleached
        else 
            bleaching_time_index(i) = NaN;%length(t_rat);
        end
    end
    bti = sort(bleaching_time_index);
    index = find(isnan(bti)==1,1,'first');
    bti(index:end)=[];
    status = ones(N,1);

    N_t = N*ones(length(t_rat),1);

    for i=1:length(bti)
       N_t(bti(i):end) = N_t(bti(i):end)-1;
    end
    
    % Plot these individual tracks
    
    I_t = N_t*I_1;
    
    % ADD NOISE?
    if add_noise==1
        noise = normrnd((I_1).*N_t/pn,I_t.*sqrt(N_t)/pn);
        % Add noise
        I_t = I_t + noise;
    end
    if add_bg_noise==1
            bg_noise = normrnd(bg,bg/10,length(I_t),1);
        I_t = I_t +bg_noise;
    end

   subplot(ceil(sqrt(T)),ceil(sqrt(T)),t)
 plot(t_rat,I_t)
    xlabel('t/\tau ratio')
    ylabel('Intensity')
    
    
end

%  plot(t_rat,I_t)
%     xlabel('t/\tau ratio')
%     ylabel('Intensity')


% Now, following SuppMeth4 in Marks Nature paper, get the pairwise differences
% for all tracks recorded for each seperate cell.
M = round(T_max/T_res);

diff_array_all = zeros(T*(M*(M-1)/2),1);
for t=1:T
    disp(['Processing track' num2str(t)])
    diff_array = zeros(M*(M-1)/2,1);
    count = 1;
    for i=1:length(I_t)
        for j=i+1:length(I_t)

            diff_array(count) = I_t(i)-I_t(j);
            count=count+1;
            
        end
    end
%     diff_array_all((M*(M-1)/2)*(t))
%     size(diff_array_all((M*(M-1)/2)*(t-1)+1:(M*(M-1)/2)*t));
    diff_array_all((M*(M-1)/2)*(t-1)+1:(M*(M-1)/2)*t)=diff_array;
%     diff_array_all((M*(M-1)/2)*(t))
    
    
%     %diff_array = diff_array./((N*(N-1)/2));
% %     subplot(3,1,2)
%     subplot(ceil(sqrt(T)),ceil(sqrt(T)),t)
%     [Nhist,X] = hist(diff_array,3000);
%     Nhist = Nhist./((N*(N-1)/2));
%     bar(X,Nhist,'hist')
%     xlabel('Pairwise intensity difference')
%     ylabel('Normalised Number of observations')
   
end

% Combine these data
    %diff_array = diff_array./((N*(N-1)/2));
%     subplot(3,1,2)
figure
    
    [Nhist,X] = hist(diff_array_all,300);
    Nhist = Nhist./((S*(S-1)/2));
    % remove all the negative elements
    [I] = find(X<0,1,'last');
    X(1:I) = [];
    Nhist(1:I) = [];
    bar(X,Nhist,'hist')
    xlabel('Pairwise intensity difference')
    ylabel('Normalised Number of observations')

%CHECK PERFORMED: THE TOTAL NUMBER OF POINTS IN DIFF_ARRAY_ALL is
%T*(M*M-1)/2) where M is the number of data points in each track with is
%thus round(T_max/T_res)

    bin_difference = diff(X);
    Fs = 1/bin_difference(1);
    L = 300;
    % S is the number of samples
    NFFT = 2^nextpow2(L);
    F = Fs/2*linspace(0,1,NFFT/2-1);
    Y = fft(Nhist,NFFT)/L;
    plot(F,abs(Y(1:NFFT/2-1)))
    Ys = abs(Y(1:NFFT/2-1));
    [I,J] = max(Ys);
    [pks,locs] = findpeaks(Ys);
    hold on
    stem(F(locs),Ys(locs),'r')
    % add in the 'real' frequency
    stem(1/I_1,0.05,'g')
    disp(['Intensity of the unitary step is' num2str(1./F(J))])





 figure
% 
pairwise_frequency = abs(ifft(Nhist));
plot(X,pairwise_frequency)
xlabel('Pairwise intensity Difference')
ylabel('FFT of frequency of observations of pairwise intensity difference')

[I,J] = max(pairwise_frequency);
disp(['Intensity of the unitary step is' num2str(X(J))])

% DIAGNOSIS OF CODE - try this with the bleaching time index




