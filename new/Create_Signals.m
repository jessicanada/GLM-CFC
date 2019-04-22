modulation_level = 1;

dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
%N  = 20/dt+4000;                           % # steps to simulate.
N  = 400/dt+4000; N = N/2;
t = (1:N)*dt;                               % Duration of simulation [s].

VpinkTest = make_pink_noise(1,N,dt);         
VpinkTest = VpinkTest - mean(VpinkTest);

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
VLOW = filtfilt(filtwts,1,Vpink);            % Define low freq band activity.
VLOW = VLOW(2001:end-2000);

% Make the data.
Vpink = make_pink_noise(1,N,dt);            % First, make pink noise.
Vpink = Vpink - mean(Vpink);                % ... with zero mean.

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
Vhi = Vhi(2001:end-2000);
t = t(2001:end-2000);

% Find peaks AND TROUGHS of the low freq activity.
[pks, ipks] = findpeaks((VLOW));          % NOTE use of "abs" here.

s = zeros(size(Vhi));                               % Define empty modulation envelope.
for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        s(ipks(i0)-10:ipks(i0)+10) = hann(21);      % ... bump up modulation envelope with Hann window.
    end
end
VHI = Vhi.*(1+modulation_level*s);                  % Old Method

%Ahi = (1+modulation_level*s)';                  % Define Ahi.
%VHI = (0.01* Ahi .* cos(angle(hilbert(Vhi'))))';% ... and use PhiHi to get Vhi.