function [XX,P,Vlo,Vhi,t] = simfun(pac_mod,aac_mod,sim_method,pval,ci,varargin)
%INPUTS:
% pac_mod:          Intensity of PAC (I_PAC in paper)
% aac_mod:          Intensity of AAC (I_AAC in paper)
% sim_method:       'GLM' creates coupled signals with a GLM approach,
%                   'pink' creates coupled signals by filtering pink noise
% pval:          	'theoretical' gives analytic p-values for R
%                   'empirical' gives bootstrapped p-values for R
% ci:               'ci' gives confidence intervals for R
%                   'none' gives no confidence intervals (faster)
% varargin:         optionally, include the parameter q indicating which quantiles
%                   of AmpLo you'd like to fit over
%
%OUTPUTS:
% XX.rpac:      R_PAC value, confidence intervals XX.rPAC_CI
% XX.raac:      R_AAC value, confidence intervals XX.rAAC_CI
% XX.rcfc:      R_CFC value, confidence intervals XX.rCFC_CI
% XX.null:      3D surface for null model in Phi_low, A_low, A_high space
% XX.PAC:       3D surface for PAC model in Phi_low, A_low, A_high space
% XX.AAC:       3D surface for AAC model in Phi_low, A_low, A_high space
% XX.CFC:       3D surface for CFC model in Phi_low, A_low, A_high space
% P.rpac:       p-value for RPAC statistic
% P.raac:       p-value for RAAC statistic
% P.rcfc:       p-value for RCFC statistic
% Vlo:          simulated low-frequency signal
% Vhi:          simulated high-frequency signal
% t:            time

dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
N  = 20/dt+4000;                            % # steps to simulate, making the duration 20s

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

% Make the data.
Vpink = make_pink_noise(1,N,dt);            % Regenerate the pink noise.
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
Vlo = Vlo(2001:end-2000);
Vhi = Vhi(2001:end-2000);
t   = (1:length(Vlo))*dt;
N   = length(Vlo);

% Find peaks of the low freq activity.
[~, ipks] = findpeaks(Vlo);
AmpLo = abs(hilbert(Vlo));

s = zeros(size(Vhi));                               % Define empty modulation envelope.
for i0=1:length(ipks)                               % At every low freq peak,
    if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
        s(ipks(i0)-10:ipks(i0)+10) = hann(21); % Scaled Modulation
    end
end
s = s/max(s);                                       % Normalize so it falls between 0 and 1

if exist('sim_method','var') && strcmp(sim_method, 'GLM')
    Ahi = (1+pac_mod*s)';                           % Define Ahi.
    Vhi = (0.01* Ahi .* cos(angle(hilbert(Vhi'))))';% ... and use PhiHi to get Vhi.
elseif exist('sim_method','var') && strcmp(sim_method, 'pink')
    Vhi = Vhi.*(1+pac_mod*s);                       % Modulate high freq activity by modulation envelope.
else
    return
end

Vhi = Vhi.*(1+aac_mod*AmpLo/max(AmpLo));

Vpink2 = make_pink_noise(1,N,dt);                   % Create additional pink noise signal
noise_level = 0.01;
V1 = Vlo+Vhi+noise_level*Vpink2;                 	% Create final "observed" new signal for analysis.

%Filter the new signal into low and high frequencies

%Filter into low freq band
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

if isempty(varargin)
    [XX,P] = glmfun(Vlo, Vhi, pval,ci);
else
    q = varargin{1};
    [XX,P] = glmfun(Vlo, Vhi, pval,ci,q);
end

end

function [x1new] = make_pink_noise(alpha,L,dt)

  %alpha=0.33;

  x1 = randn(L,1);
  xf1 = fft(x1);
  A = abs(xf1);
  phase = angle(xf1);

  df = 1.0 / (dt*length(x1));
  faxis = (0:length(x1)/2)*df;
  faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
  oneOverf = 1.0 ./ faxis.^alpha;
  oneOverf(1)=0.0;

  Anew = A.*oneOverf';
  xf1new = Anew .* exp(i*phase);
  x1new = real(ifft(xf1new))';
  
end