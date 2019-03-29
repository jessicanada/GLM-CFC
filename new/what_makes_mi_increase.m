% 
% %% What makes MI increase/decrease?
% N = 100;
% mi_change = zeros(1,N);
% rpac_change = zeros(1,N);
% for i = 1:N
% % to simulate: positive coupling at pi/2, negative at -pi/2? or positive at
% % pi negative at 0, whatever
% dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
% N  = 20/dt+4000;
% 
% % Filter into low freq band.
% locutoff = 4;                               % Low freq passband = [4,7] Hz.
% hicutoff = 7;
% filtorder = 3*fix(Fs/locutoff);
% MINFREQ = 0;
% trans          = 0.15;                      % fractional width of transition zones
% f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
% m=[0       0                      1            1            0                      0];
% filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients
% 
% % Filter into high freq band.
% locutoff = 100;                             % High freq passband = [100, 140] Hz.
% hicutoff = 140;
% filtorder = 10*fix(Fs/locutoff);
% MINFREQ = 0;
% trans          = 0.15;                      % fractional width of transition zones
% f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
% m=[0       0                      1            1            0                      0];
% filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients
%                           % # steps to simulate, making the duration 20s
% 
% % Make the data.
% Vpink = make_pink_noise(1,N,dt);            % First, make pink noise.
% Vpink = Vpink - mean(Vpink);                % ... with zero mean.
% Vlo = filtfilt(filtwts_lo,1,Vpink);            % Define low freq band activity.
% 
% % Make the data.
% Vpink = make_pink_noise(1,N,dt);            % Regenerate the pink noise.
% Vpink = Vpink - mean(Vpink);                % ... with zero mean.
% Vhi = filtfilt(filtwts_hi,1,Vpink);            % Define high freq band activity.
% 
% % Drop the edges of filtered data to avoid filter artifacts.
% Vlo = Vlo(2001:end-2000);
% Vhi = Vhi(2001:end-2000);
% t   = (1:length(Vlo))*dt;
% N   = length(Vlo);
% 
% % Find peaks of the low freq activity.
% 
% AmpLo = abs(hilbert(Vlo));
% 
% %monophasic coupling
% [~, ipks] = findpeaks((Vlo));
% s_mono = zeros(1,length(Vhi));                               % Define empty modulation envelope.
% 
% for i0=1:length(ipks)                               % At every low freq peak,
%     if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
%         s_mono(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
%     end
% end
% 
% s_mono = s_mono/max(s_mono);
% 
% %weird coupling in POST
% [~, ipks] = findpeaks(abs(Vlo));
% s_post = zeros(1,length(Vhi));                               % Define empty modulation envelope.
% 
% for i0=1:length(ipks)                               % At every low freq peak,
%     if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
%         if Vlo(ipks(i0))>0
%             s_post(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
%         else
%             s_post(ipks(i0)-10:ipks(i0)+10) = -hann(21); % if it's a trough
%         end
%     end
% end
% 
% s_post = s_post/max(s_post);
% 
% %normal coupling in PRE
% [~, ipks] = findpeaks(abs(Vlo));
% s_pre = zeros(1,length(Vhi));
% for i0=1:length(ipks)                               % At every low freq peak,
%     if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
%         if Vlo(ipks(i0))>0
%             s_pre(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
%         else
%             s_pre(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a trough
%         end
%     end
% end
% 
% %s = [s_pre(1:length(s_pre)/2),s_post(length(s_pre)/2+1:end)];
% s = s_mono;
% 
% aac_mod = [0*ones(length(Vhi)/2,1);5*ones(length(Vhi)/2,1)]';
% pac_mod = [0*ones(length(Vhi)/2,1);0*ones(length(Vhi)/2,1)]'; %increase PAC in post
% Vhi = Vhi.*(1+pac_mod.*s);                       
% Vhi = Vhi.*(1+aac_mod.*AmpLo/max(AmpLo));
% 
% Vpink2 = make_pink_noise(1,N,dt);                   
% noise_level = 0.01;
% 
% Vlo = Vlo.*[1*ones(length(Vhi)/2,1);1*ones(length(Vhi)/2,1)]'; %increase low frequency amplitude in post
% Vhi = Vhi.*[1*ones(length(Vhi)/2,1);1*ones(length(Vhi)/2,1)]'; %increase high frequency amplitude in post
% V1 = Vlo+Vhi+noise_level*Vpink2;                 	
% 
% %Filter into low freq band
% Vlo = filtfilt(filtwts_lo,1,V1);            % Define low freq band activity.
% 
% % Filter into high freq band.
% Vhi = filtfilt(filtwts_hi,1,V1);            % Define high freq band activity.
% 
% [XX1] = glmfun(Vlo(1:length(Vhi)/2), Vhi(1:length(Vhi)/2), 'none','none',.05);
% [MI1] = modulation_index(Vlo(1:length(Vhi)/2),Vhi(1:length(Vhi)/2),'none');
% 
% [XX2] = glmfun(Vlo(length(Vhi)/2:end), Vhi(length(Vhi)/2:end), 'none','none',.05);
% [MI2] = modulation_index(Vlo(length(Vhi)/2:end),Vhi(length(Vhi)/2:end),'none');
% 
% if MI1 < MI2
%     mi_change(i) = 1;
% else
%     mi_change(i) = 0;
% end
%     
% if XX1.rpac_new < XX2.rpac_new
%     rpac_change(i) = 1;
% else
%     rpac_change(i) = 0;
% end
% 
% if XX1.rpac_new_mak < XX2.rpac_new_mak
%     rpac_change_mak(i) = 1;
% else
%     rpac_change_mak(i) = 0;
% end
% 
% end
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp(['Mean MI Change: ',num2str(mean(mi_change))])
% %disp(['MI CI: [',num2str(mean(mi_change)-1.96*std(mi_change)),', ',num2str(mean(mi_change)+1.96*std(mi_change)),']'])
% disp(['Mean RPAC Change: ',num2str(mean(rpac_change))])
% 
% %%
% figure(2)
% subplot(1,2,1)
% modulation_index(Vlo(1:length(Vhi)/2),Vhi(1:length(Vhi)/2),'plot');
% title('pre')
% ylim([0,.09])
% subplot(1,2,2)
% modulation_index(Vlo(length(Vhi)/2:end),Vhi(length(Vhi)/2:end),'plot');
% title('post')
% ylim([0,.09])
% % modulation index is decreasing
% %%
% figure(1)
% subplot(1,2,1)
% [mi_pre,a_mean_pre] = modulation_index(Vlo_pre,Vhi_pre,'plot');
% title('pre')
% ylim([.05,.06])
% subplot(1,2,2)
% [mi_post,a_mean_post] = modulation_index(Vlo_post,Vhi_post,'plot');
% title('post')
% ylim([.05,.06])
% % modulation index is decreasing
% %%
% %idea: look at chunks of pre vs. chunks of post?
% sum(abs(a_mean_pre-1/18))
% sum(abs(a_mean_post-1/18))
% 

%% What makes MI increase/decrease?
N = 100;
RPAC1 = zeros(1,N); p_RPAC1 = zeros(1,N);
MI1 = zeros(1,N); p_MI1 = zeros(1,N);
RPAC2 = zeros(1,N); p_RPAC2 = zeros(1,N);
MI2 = zeros(1,N); p_MI2 = zeros(1,N);
for i = 1:N
    i
% to simulate: positive coupling at pi/2, negative at -pi/2? or positive at
% pi negative at 0, whatever
dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
N  = 20/dt+4000;

% Filter into low freq band.
locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 7;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients

% Filter into high freq band.
locutoff = 100;                             % High freq passband = [100, 140] Hz.
hicutoff = 140;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients
                          % # steps to simulate, making the duration 20s

% Make the data.
Vpink = make_pink_noise(1,N,dt);            % First, make pink noise.
Vpink = Vpink - mean(Vpink);                % ... with zero mean.
Vlo = filtfilt(filtwts_lo,1,Vpink);            % Define low freq band activity.

% Make the data.
Vpink = make_pink_noise(1,N,dt);            % Regenerate the pink noise.
Vpink = Vpink - mean(Vpink);                % ... with zero mean.
Vhi = filtfilt(filtwts_hi,1,Vpink);            % Define high freq band activity.

% Drop the edges of filtered data to avoid filter artifacts.
Vlo = Vlo(2001:end-2000);
Vhi = Vhi(2001:end-2000);
t   = (1:length(Vlo))*dt;
N   = length(Vlo);

% Find peaks of the low freq activity.

AmpLo = abs(hilbert(Vlo));

%monophasic coupling
[~, ipks] = findpeaks((Vlo));
s_mono = zeros(1,length(Vhi));                               % Define empty modulation envelope.

for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        s_mono(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
    end
end

s_mono = s_mono/max(s_mono);
s = s_mono;

aac_mod = 5*ones(length(Vhi),1)';
pac_mod = 0*ones(length(Vhi),1)'; %increase PAC in post
Vhi = Vhi.*(1+pac_mod.*s);                       
Vhi1 = Vhi.*(1+aac_mod.*AmpLo/max(AmpLo)); %with AAC
Vhi2 = Vhi; %without AAC

Vpink2 = make_pink_noise(1,N,dt);                   
noise_level = 0.01;

%Vlo = Vlo.*[1*ones(length(Vhi)/2,1);1*ones(length(Vhi)/2,1)]'; %increase low frequency amplitude in post
%Vhi = Vhi.*[1*ones(length(Vhi)/2,1);1*ones(length(Vhi)/2,1)]'; %increase high frequency amplitude in post
V1 = Vlo+Vhi1+noise_level*Vpink2; 
V2 = Vlo+Vhi2+noise_level*Vpink2;

%Filter into low freq band
Vlo1 = filtfilt(filtwts_lo,1,V1);            % Define low freq band activity.
% Filter into high freq band.
Vhi1 = filtfilt(filtwts_hi,1,V1);            % Define high freq band activity.

%Filter into low freq band
Vlo2 = filtfilt(filtwts_lo,1,V2);            % Define low freq band activity.
% Filter into high freq band.
Vhi2 = filtfilt(filtwts_hi,1,V2);            % Define high freq band activity.

%with AAC
[XX1,P1] = glmfun(Vlo1, Vhi1, 'empirical','none',.05);
[MI1,MP1] = modulation_index(Vlo1,Vhi1,'pvals');
RPAC1(i) = XX1.rpac_new; p_RPAC1(i) = P1.rpac_new;
MI1(i) = MI1; p_MI1(i) = MP1;

%without AAC
[XX2,P2] = glmfun(Vlo2, Vhi2, 'empirical','none',.05);
[MI2,MP2] = modulation_index(Vlo2,Vhi2,'pvals');
RPAC2(i) = XX2.rpac_new; p_RPAC2(i) = P2.rpac_new;
MI2(i) = MI2; p_MI2(i) = MP2;
end

save('R_MI_Comparison_Increase_AAC','RPAC1','RPAC2','p_RPAC1','p_RPAC2','MI1','MI2','p_MI1','p_MI2')