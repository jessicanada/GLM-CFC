addpath('Chaotic Systems Toolbox')

aac = 1:11;
RPAC = zeros(length(aac),1000); MI = zeros(length(aac),1000);
for j = 1:1000
    j
    pac_mod = 1;
    dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
    N  = 400/dt+4000; N = N/2;                          % # steps to simulate, making the duration 20s

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
    filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients

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

    Vhi = filtfilt(filtwts_hi,1,Vpink);            % Define high freq band activity.

    % Drop the edges of filtered data to avoid filter artifacts.
    Vlo = Vlo(2001:end-2000); VLO = Vlo;
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


    VHI = Vhi.*(1+pac_mod*s);                       % Modulate high freq activity by modulation envelope.
    Vpink2 = make_pink_noise(1,N,dt);                   % Create additional pink noise signal
    for i = 1:length(aac)
        Vhi1 = VHI.*(1+aac(i)*AmpLo/max(AmpLo));

        noise_level = 0.01;
        V1 = VLO+Vhi1+noise_level*Vpink2;                 	% Create final "observed" new signal for analysis.

        Vlo = filtfilt(filtwts_lo,1,V1);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V1);            % Define high freq band activity.

        [XX] = glmfun(Vlo,Vhi,'none','none','none',.05);
        RPAC(i,j) = XX.rpac_new;
        [mi] = modulation_index(Vlo,Vhi,'none');
        MI(i,j) = mi;
    end
end
strname = ['Increase_AAC_Results_200_Seconds_'];
save(strname,'RPAC','MI')