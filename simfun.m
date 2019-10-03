function [XX,P,Vlo,Vhi,t] = simfun(pac_mod,aac_mod,sim_method,pval,ci,AIC,varargin)
%INPUTS:
% pac_mod:          Intensity of PAC (I_PAC in paper)
% aac_mod:          Intensity of AAC (I_AAC in paper)
% sim_method:       'GLM' creates coupled signals with a GLM approach,
%                   'pink' creates coupled signals by filtering pink noise
% pval:          	'theoretical' gives analytic p-values for R
%                   'empirical' gives bootstrapped p-values for R
% ci:               'ci' gives confidence intervals for R
%                   'none' gives no confidence intervals (faster)
% AIC:              'AIC' computes number of control points for spline
%                   phase via AIC minimization
% varargin:         optionally, include the parameter q indicating which quantiles
%                   of AmpLo you'd like to fit over
%
%OUTPUTS:
% XX.rpac:                R_PAC value, confidence intervals XX.rPAC_CI
% XX.raac:                R_AAC value, confidence intervals XX.rAAC_CI
% XX.Phi_low:             3D surface for Phi_low model in Phi_low, A_low, A_high space
% XX.A_low:               3D surface for A_low model in Phi_low, A_low, A_high space
% XX.Phi_low_A_low:       3D surface for Phi_low,A_low model in Phi_low, A_low, A_high space
% P.rpac:                 p-value for RPAC statistic
% P.raac:                 p-value for RAAC statistic
% Vlo:                    simulated low-frequency signal
% Vhi:                    simulated high-frequency signal
% t:                      time

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

if exist('sim_method','var')
    
    switch sim_method
        
        case 'GLM'
            Ahi = (1+pac_mod*s)';                           % Define Ahi.
            Vhi = (0.01* Ahi .* cos(angle(hilbert(Vhi'))))';% ... and use PhiHi to get Vhi.
            Vhi = Vhi.*(1+aac_mod*AmpLo/max(AmpLo));
            
        case 'pink'
            Vhi = Vhi.*(1+pac_mod*s);                       % Modulate high freq activity by modulation envelope.
            Vhi = Vhi.*(1+aac_mod*AmpLo/max(AmpLo));
      
        case 'spiking'                                      % Add a spiking model.
            N     = 20/dt;                                  % # steps to simulate, making the duration 20s
            t     = (1:N)*dt;                               % Time axis.
            Alo   = 1+(sin(2*pi*t*0.1)+1)/2;                % Slow modulation of low frequency envelope.
            Philo = pi*sawtooth(2*pi*t*4);                  % Low frequency phase is periodic (4 Hz).
            Vlo   = Alo.*cos(Philo);

            Philo_star = pi + Alo*pi;                       % Target phase depends on low frequency envelope.
            sigma = 0.01;
            lambda = 1/sqrt(2*pi*sigma) * exp( -(1+sawtooth(Philo - Philo_star,1/2)).^2 / (2*sigma^2) );
            lambda = (0.001+0.3*lambda/max(lambda));
            Vhi = 1*binornd(1,lambda);                      % When low freq phase is near target phase, produce a spike.
            
            Vlo =(Vlo)+0.1*randn(size(Vlo));                % Define Vlo and Vhi directly for this sim.
            Vhi =(Vhi)+0.1*randn(size(Vlo));
    end
    
else
    return
end

if exist('sim_method','var') && ~strcmp(sim_method, 'spiking')

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

end

if isempty(varargin)
    [XX,P] = glmfun(Vlo, Vhi, pval,ci,AIC);
else
    q = varargin{1};
    [XX,P] = glmfun(Vlo, Vhi, pval,ci,AIC,q);
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