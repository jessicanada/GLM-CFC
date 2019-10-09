%% Simulate data for pre/post comparison
function [Vlo_pre,Vhi_pre,Vlo_post,Vhi_post] = simulate_prepost_increase_alow(pac_mod_1,pac_mod_2,aac_mod_1,aac_mod_2)
    %pac_mod_1 = 0;
    %pac_mod_2 = 3;
    %aac_mod_1 = 1;
    %aac_mod_2 = 1;

    dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
    N  = 40/dt+4000;                            % # steps to simulate, making the duration 20s

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
    filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients
    Vhi = filtfilt(filtwts_hi,1,Vpink);            % Define high freq band activity.

    % Drop the edges of filtered data to avoid filter artifacts.
    Vlo = Vlo(2001:end-2000);
    Vlo = Vlo.*[ones(1,length(Vlo)/2),2*ones(1,length(Vlo)/2)];
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

    PAC = [pac_mod_1*ones(1,length(s)/2),pac_mod_2*ones(1,length(s)/2)];
    AAC = [aac_mod_1*ones(1,length(s)/2),aac_mod_2*ones(1,length(s)/2)];

    Vhi = Vhi.*(1+PAC.*s);                       % Modulate high freq activity by modulation envelope.
    Vhi = Vhi.*(1+AAC.*AmpLo/max(AmpLo));

    Vpink2 = make_pink_noise(1,N,dt);                   % Create additional pink noise signal
    noise_level = 0.01;
    V1 = Vlo+Vhi+noise_level*Vpink2;                 	% Create final "observed" new signal for analysis.

    Vlo = filtfilt(filtwts_lo,1,V1);            % Define low freq band activity.
    Vhi = filtfilt(filtwts_hi,1,V1);            % Define high freq band activity.

    Vlo_pre = Vlo(1:length(Vlo)/2)'; 
    Vhi_pre = Vhi(1:length(Vlo)/2)'; 

    Vlo_post = Vlo(length(Vlo)/2+1:end)';
    Vhi_post = Vhi(length(Vlo)/2+1:end)';

%%[XX,P,I] = glmfun_with_indicator_update(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,'none','none',.05);
%%plot([Vlo_pre;Vlo_post]); hold on; plot([Vhi_pre;Vhi_post])