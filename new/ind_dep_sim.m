function [p,xx,P,XX,mi,p_mi,MI,P_MI] = ind_dep_sim(threshold, modulation_level)

%% Create Vhi1 (dependent on AmpLo) and Vhi2 (not dependent on AmpLo)

dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
N  = 2*(30000+4000);%30000+4000;            % # steps to simulate.
N = round(N/6.8);                           % 20s duration
t = (1:N)*dt;                               % Duration of simulation [s].

% Make the data.
Vpink = make_pink_noise(1,N,dt);            % First, make pink noise.
Vpink = Vpink - mean(Vpink);                % ... with zero mean.

% Filter into low freq band.
locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 7;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients
Vlo = filtfilt(filtwts_lo,1,Vpink);            % Define low freq band activity.

% Filter into high freq band.
locutoff = 100;                             % High freq passband = [100, 140] Hz.
hicutoff = 140;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients
Vhi = filtfilt(filtwts_hi,1,Vpink);            % Define high freq band activity.

% Drop the edges of filtered data to avoid filter artifacts.
Vlo = Vlo(2001:end-2000);
Vhi = Vhi(2001:end-2000);
t   = (1:length(Vlo))*dt;

Alo = abs(hilbert(Vlo));

% Find peaks of the low freq activity.
[~, ipks] = findpeaks(Vlo);

count = 0;
s = zeros(size(Vhi));                               % Define empty modulation envelope.
for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        val = Alo(ipks(i0));
        if val>=quantile(Alo(ipks),threshold) %If the amplitude is in the top half of the distribution
            scalar = 1;               %Make the modulation large
            count = count+1;
        else
            scalar = 0;
        end
        s(ipks(i0)-10:ipks(i0)+10) = scalar*hann(21); %Scaled Modulation
    end
end
s = s/max(s); %normalize so it falls between 0 and 1
Vhi1 = Vhi.*(1+modulation_level*s);

ind = randi(length(ipks),count,1); %Choose the same number of large CFC events as before, sampled from ipks values

s2 = zeros(size(Vhi));                               % Define empty modulation envelope.
for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        if any(i0==ind) %If the amplitude is in the top half of the distribution
            scalar = 1;               %Make the modulation large
            count = count+1;
        else
            scalar = 0; %Otherwise, make the modulation small (cube of small value<1)
        end
        s2(ipks(i0)-10:ipks(i0)+10) = scalar*hann(21); %Scaled Modulation
    end
end
s2 = s2/max(s2); %normalize so it falls between 0 and 1
Vhi2 = Vhi.*(1+modulation_level*s2);

% Create VLO1/VHI1 (Dependence on AmpLo) and VLO2/VHI2 (No Dependence on AmpLo)

nCtlPts = 10;

Vpink2 = make_pink_noise(1,N,dt);
Vpink2 = Vpink2(2001:end-2000);
Vlo = Vlo+.01*Vpink2; %add noise to low-frequency signal

V1 = Vlo+Vhi1;
V2 = Vlo+Vhi2;

%recalculate Vlo, Vhi from composite signal
    
    VLO1 = filtfilt(filtwts_lo,1,V1); %dependent
    VHI1 = filtfilt(filtwts_hi,1,V1);
    
    VLO2 = filtfilt(filtwts_lo,1,V2); %independent
    VHI2 = filtfilt(filtwts_hi,1,V2);

[xx,p] = glmfun(VLO1, VHI1,'empirical','none',.05); %dependent
 [mi,p_mi] = modulation_index(VLO1,VHI1,'pvals');
[XX,P] = glmfun(VLO2, VHI2, 'empirical','none',.05); %independent
 [MI,P_MI] = modulation_index(VLO2,VHI2,'pvals');
