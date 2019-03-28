%% R, MI values (day by electrode) for single electrode data

load('Rat_Data_Pre')
RPAC_pre = RPAC; RPAC_new_pre = RPAC_new; RAAC_pre = RAAC; RAAC_new_pre = RAAC_new; RCFC_pre = RCFC; MI_pre = MI;
load('Rat_Data_Post')
RPAC_diff = sign(RPAC-RPAC_pre); RAAC_diff = sign(RAAC-RAAC_pre); RPAC_new_diff = sign(RPAC_new-RPAC_new_pre); RAAC_new_diff = sign(RAAC_new-RAAC_new_pre);
RCFC_diff = sign(RCFC-RCFC_pre); MI_diff = sign(MI-MI_pre);

%% extract data from day 2, electrode 1
clear;
load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;


dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 
locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients

locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients


k = 2;%day

j = 1; %il trial

il = ildata{k,1}; %pre
il = il(:,j);
il = decimate(il,10);
il = decimate(il,3);

Vlo_pre = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
Vhi_pre = filtfilt(filtwts_hi,1,il);

il = ildata{k,2}; %post
il = il(:,j);
il = decimate(il,10);
il = decimate(il,3);
Vlo_post = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
Vhi_post = filtfilt(filtwts_hi,1,il);

% am = amdata{k,1}; %pre
% am = am(:,j-4);
% am = decimate(am,10);
% am = decimate(am,3);
% 
% Vlo_pre = filtfilt(filtwts_lo,1,am);            % Define low freq band activity.
% Vhi_pre = filtfilt(filtwts_hi,1,am);
% 
% am = amdata{k,2}; %post
% am = am(:,j-4);
% am = decimate(am,10);
% am = decimate(am,3);
% 
% Vlo_post = filtfilt(filtwts_lo,1,am);            % Define low freq band activity.
% Vhi_post = filtfilt(filtwts_hi,1,am);
   


%% calculate R, MI values
[R_pre] = glmfun(Vlo_pre',Vhi_pre','none','none',.05);
[mi_pre] = modulation_index(Vlo_pre,Vhi_pre,'none');

[R_post] = glmfun(Vlo_post',Vhi_post','none','none',.05);
[mi_post] = modulation_index(Vlo_post,Vhi_post,'none');

[R_ind,~,I] = glmfun_with_indicator(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,'none',.05);

%AAC appears to have increased
%rpac_new says PAC increased, RPAC and MI say decreased, I says increased
%%
plot([Vlo_pre;Vlo_post]); hold on; plot([Vhi_pre;Vhi_post]); plot([length(Vlo_pre),length(Vlo_pre)],[-100,100],'k','LineWidth',2)

%%
figure(1)
histogram(abs(hilbert(Vlo_pre)),'Normalization','Probability'); hold on; histogram(abs(hilbert(Vlo_post)),'Normalization','Probability');
legend('pre','post'); title('Vlo')

figure(2)
histogram(abs(hilbert(Vhi_pre)),'Normalization','Probability'); hold on; histogram(abs(hilbert(Vhi_post)),'Normalization','Probability');
legend('pre','post'); title('Vhi')

%both low and high frequency amplitudes increased

%%

%the coupling changed! 

XX1 = R_pre;
surf(XX1.ampAXIS,XX1.phi0,XX1.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX1.ampAXIS,XX1.phi0,XX1.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX1.ampAXIS,XX1.phi0,XX1.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
surf(XX1.ampAXIS,XX1.phi0,XX1.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX1.ampAXIS),max(XX1.ampAXIS)]); ylim([min(XX1.phi0),max(XX1.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('PAC','AAC','CFC')
title('pre')
set(gca,'FontSize',13)
grid off

figure(2)
XX2 = R_post;
surf(XX2.ampAXIS,XX2.phi0,XX2.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX2.ampAXIS,XX2.phi0,XX2.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX2.ampAXIS,XX2.phi0,XX2.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
surf(XX2.ampAXIS,XX2.phi0,XX2.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX2.ampAXIS),max(XX2.ampAXIS)]); ylim([min(XX2.phi0),max(XX2.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('PAC','AAC','CFC')
title('post')
set(gca,'FontSize',13)
grid off

%% overall
[R] = glmfun([Vlo_pre',Vlo_post'],[Vhi_pre',Vhi_post'],'none','none',.05);
XX1 = R;
surf(XX1.ampAXIS,XX1.phi0,XX1.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX1.ampAXIS,XX1.phi0,XX1.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX1.ampAXIS,XX1.phi0,XX1.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
surf(XX1.ampAXIS,XX1.phi0,XX1.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX1.ampAXIS),max(XX1.ampAXIS)]); ylim([min(XX1.phi0),max(XX1.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('PAC','AAC','CFC')
title('pre')
set(gca,'FontSize',13)
grid off
%%
figure(1)
subplot(1,2,1)
modulation_index(Vlo_pre,Vhi_pre,'plot')
title('pre')
subplot(1,2,2)
modulation_index(Vlo_post,Vhi_post,'plot')
title('post')
%% try and recreate this in simulation
N = 100;
MI_change = zeros(1,N);
RPAC_change = zeros(1,N);
RPAC_new_change = zeros(1,N);

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

    %weird coupling in POST
    [~, ipks] = findpeaks(abs(Vlo));
    s_post = zeros(1,length(Vhi));                               % Define empty modulation envelope.

    for i0=1:length(ipks)                               % At every low freq peak,
        if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
            if Vlo(ipks(i0))>0
                s_post(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
            else
                s_post(ipks(i0)-10:ipks(i0)+10) = -hann(21); % if it's a trough
            end
        end
    end

    s_post = s_post/max(s_post);

    %normal coupling in PRE
    [~, ipks] = findpeaks(abs(Vlo));
    s_pre = zeros(1,length(Vhi));
    for i0=1:length(ipks)                               % At every low freq peak,
        if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
            if Vlo(ipks(i0))>0
                s_pre(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a peak
            else
                s_pre(ipks(i0)-10:ipks(i0)+10) = hann(21); % if it's a trough
            end
        end
    end

    s = [s_pre(1:length(s_pre)/2),s_post(length(s_pre)/2+1:end)];
    aac_mod = [2*ones(length(Vhi)/2,1);3*ones(length(Vhi)/2,1)]';
    pac_mod = [.4*ones(length(Vhi)/2,1);.5*ones(length(Vhi)/2,1)]'; %increase PAC in post
    Vhi = Vhi.*(1+pac_mod.*s);                       
    Vhi = Vhi.*(1+aac_mod.*AmpLo/max(AmpLo));

    Vpink2 = make_pink_noise(1,N,dt);                   
    noise_level = 0.01;

    Vlo = Vlo.*[1*ones(length(Vhi)/2,1);5*ones(length(Vhi)/2,1)]'; %increase low frequency amplitude in post
    Vhi = Vhi.*[1*ones(length(Vhi)/2,1);5*ones(length(Vhi)/2,1)]'; %increase high frequency amplitude in post
    V1 = Vlo+Vhi+noise_level*Vpink2;                 	

    %Filter into low freq band
    Vlo = filtfilt(filtwts_lo,1,V1);            % Define low freq band activity.

    % Filter into high freq band.
    Vhi = filtfilt(filtwts_hi,1,V1);            % Define high freq band activity.

    [XX1] = glmfun(Vlo(1:length(Vhi)/2), Vhi(1:length(Vhi)/2), 'none','none',.05);
    [MI1] = modulation_index(Vlo(1:length(Vhi)/2),Vhi(1:length(Vhi)/2),'none');

    [XX2] = glmfun(Vlo(length(Vhi)/2:end), Vhi(length(Vhi)/2:end), 'none','none',.05);
    [MI2] = modulation_index(Vlo(length(Vhi)/2:end),Vhi(length(Vhi)/2:end),'none');

%     if MI1 < MI2
%         MI_change(i) = 1;
%     else
%         MI_change(i) = -1;
%     end
% 
%     if XX1.rpac_new < XX2.rpac_new
%         RPAC_new_change(i) = 1;
%     else
%         RPAC_new_change(i) = -1;
%     end
% 
%     if XX1.rpac < XX2.rpac
%         RPAC_change(i) = 1;
%     else
%         RPAC_change(i) = -1;
%     end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


if MI1 < MI2
    disp('MI Increased')
else
    disp('MI Decreased')
end
    
if XX1.rpac_new < XX2.rpac_new
    disp('RPAC_new Increased')
else
    disp('RPAC_new Decreased')
end

if XX1.rpac < XX2.rpac
    disp('RPAC Increased')
else
    disp('RPAC Decreased')
end
  
%%
figure(1)
subplot(2,1,1)
surf(XX1.ampAXIS,XX1.phi0,XX1.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX1.ampAXIS,XX1.phi0,XX1.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX1.ampAXIS,XX1.phi0,XX1.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
surf(XX1.ampAXIS,XX1.phi0,XX1.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX1.ampAXIS),max(XX1.ampAXIS)]); ylim([min(XX1.phi0),max(XX1.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
%legend('PAC','AAC','null','CFC')
title('pre')
set(gca,'FontSize',13)
grid off

subplot(2,1,2)
surf(XX2.ampAXIS,XX2.phi0,XX2.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX2.ampAXIS,XX2.phi0,XX2.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX2.ampAXIS,XX2.phi0,XX2.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
surf(XX2.ampAXIS,XX2.phi0,XX2.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX2.ampAXIS),max(XX2.ampAXIS)]); ylim([min(XX2.phi0),max(XX2.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
%legend('PAC','AAC','null','CFC')
title('post')
set(gca,'FontSize',13)
grid off

%% Code to analyze single electrode rat data
clear;
load('/Users/jnadalin/Dropbox/Nadalin (1)/Raw_data/Good_phase_locking_N=4/data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 16 combinations
RPAC_new = zeros(4,8);
RAAC = zeros(4,8); 
RAAC_new = zeros(4,8);
RCFC = zeros(4,8); 
MI = zeros(4,8);
locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients

locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients

%post
for k = 1:4 %day
        count = 1;
        for j = 1:4 %il trial
            il = ildata{k,2}; %post
            il = il(:,j);
            il = decimate(il,10);
            il = decimate(il,3);
            
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            Vhi = filtfilt(filtwts_hi,1,il);
            
            [r] = glmfun(Vlo',Vhi','none','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            RPAC(k,count) = r.rpac;
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
        for i = 1:4 %am trial
            am = amdata{k,2}; %post
            am = am(:,i);
            am = decimate(am,10);
            am = decimate(am,3);

            % Filter into high freq band.
            Vhi = filtfilt(filtwts_hi,1,am);            % Define high freq band activity.
            Vlo = filtfilt(filtwts_lo,1,am);
            [r] = glmfun(Vlo',Vhi','none','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            RPAC(k,count) = r.rpac;
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
end
save('Rat_Data_Post','RPAC','RPAC_new','RAAC','RAAC_new','RCFC','MI')

%%
clear;
load('/Users/jnadalin/Dropbox/Nadalin (1)/Raw_data/Good_phase_locking_N=4/data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 16 combinations
RPAC_new = zeros(4,8);
RAAC = zeros(4,8); 
RAAC_new = zeros(4,8);
RCFC = zeros(4,8); 
MI = zeros(4,8);

locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients

locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients

%post
for k = 1:4 %day
        count = 1;
        for j = 1:4 %il trial
            il = ildata{k,1}; %pre
            il = il(:,j);
            il = decimate(il,10);
            il = decimate(il,3);
            
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            Vhi = filtfilt(filtwts_hi,1,il);
            
            [r] = glmfun(Vlo',Vhi','none','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            RPAC(k,count) = r.rpac;
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
        for i = 1:4 %am trial
            am = amdata{k,1}; %pre
            am = am(:,i);
            am = decimate(am,10);
            am = decimate(am,3);

            % Filter into high freq band.
            Vhi = filtfilt(filtwts_hi,1,am);            % Define high freq band activity.
            Vlo = filtfilt(filtwts_lo,1,am);
            [r] = glmfun(Vlo',Vhi','none','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            RPAC(k,count) = r.rpac;
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
end
save('Rat_Data_Pre','RPAC','RPAC_new','RAAC','RAAC_new','RCFC','MI')

