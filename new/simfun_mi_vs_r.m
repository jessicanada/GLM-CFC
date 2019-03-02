function [XX,P,Vlo,Vhi,t] = simfun_mi_vs_r(modulation_level,fixed_val)
%INPUTS:
%modulation_level: Value of modulation strength in [0,1]
%threshold: Thresholding value in [0,1], e.g. .75 means only the top 25%
%highest amplitudes are modulated
%varargin: If 'ON', computes spline for different values of ampLO

%OUTPUTS:
%R: r distance values from the null model for models 1, 2, and 3
%XX: Structure of spline matrices, each column corresponding to a different
%ampLo value
%chitest: Chi^2 test, determining the significance of the interaction term
%ampLO*phiLO in model 3

dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
N  = 2*(30000+4000);%30000+4000;            % # steps to simulate.
N = round(N/6.8); 
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
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vlo = filtfilt(filtwts,1,Vpink);            % Define low freq band activity.

% Filter into high freq band.
locutoff = 100;                             % High freq passband = [100, 140] Hz.
hicutoff = 140;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vhi = filtfilt(filtwts,1,Vpink);            % Define high freq band activity.

% Drop the edges of filtered data to avoid filter artifacts.
Vlo = Vlo(2001:end-2000);
Vhi = Vhi(2001:end-2000);
t   = (1:length(Vlo))*dt;
Vhi_old = Vhi;

maxamp = max(abs(hilbert(Vhi)));

% Find peaks of the low freq activity.
[pks, ipks] = findpeaks(Vlo);
AmpLo = abs(hilbert(Vlo));
Vlo_old = Vlo;

scalar = [];
s1 = zeros(size(Vhi));
s = zeros(size(Vhi));  
THRESH = quantile(AmpLo(ipks),.5);
% Define empty modulation envelope.
for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        val = AmpLo(ipks(i0));
        if val>=quantile(AmpLo(ipks),.5) %If the amplitude is in the top half of the distribution
            scalar = 1;               %Make the modulation large
        else
            scalar = fixed_val; %Otherwise, make the modulation small (cube of small value<1)
            %scalar = 0;
        end
        s(ipks(i0)-10:ipks(i0)+10) = scalar*hann(21); %Scaled Modulation
        s1(ipks(i0)-10:ipks(i0)+10) = hann(21);% Unscaled Modulation
    end
end
s = s/max(s); %normalize so it falls between 0 and 1
%Ahi = (1+modulation_level*s)';                  % Define Ahi.
%Vhi = (0.01* Ahi .* cos(angle(hilbert(Vhi'))))';% ... and use PhiHi to get Vhi.
Vhi = Vhi.*(1+modulation_level*s); %


Vpink2 = make_pink_noise(1,N,dt);
Vpink2 = Vpink2(2001:end-2000);
Vlo = Vlo+.01*Vpink2;

V1 = Vlo+Vhi;
%figure(2)
%plot(t,V1,t,Vlo,t,Vhi); legend('V1','Vlo','Vhi')


locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 7;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vlo = filtfilt(filtwts,1,V1);            % Define low freq band activity.

% Filter into high freq band.
locutoff = 100;                             % High freq passband = [100, 140] Hz.
hicutoff = 140;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vhi = filtfilt(filtwts,1,V1);            % Define high freq band activity.

%     dt = t(3)-t(2);
%     fnq = 1/dt/2;
%     n = 500;
%     wn = [4/fnq;7/fnq];
%     b = transpose(fir1(n,wn,'bandpass'));
%     Vlo = filtfilt(b,1,V1);
%     
%     wn = [100/fnq;140/fnq]; %85 gets rid of the curve that 70 has 70 250
%     b = transpose(fir1(n,wn,'bandpass'));
%     Vhi = filtfilt(b,1,V1);
%figure(3)
%plot(t,V1,t,Vlo,t,Vhi); legend('V1','Vlo','Vhi')

[XX,P] = glmfun(Vlo, Vhi,'empirical','none',.05);
%[XX,P] = glmfun(Vlo, Vhi,'empirical',.05); %Changed to empirical p-values (06/19/18)
end