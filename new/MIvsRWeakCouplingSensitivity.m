addpath('new')
addpath('Chaotic Systems Toolbox')

PAC_mod = 0:0.05:.5;
J = 1000;
rcfc = zeros(length(PAC_mod),J); rpac = zeros(length(PAC_mod),J);
p_rcfc = zeros(length(PAC_mod),J); p_rpac = zeros(length(PAC_mod),J);
mi = zeros(length(PAC_mod),J); p_mi = zeros(length(PAC_mod),J);
for j = 1:J
    j
    for i = 1:length(PAC_mod)
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

        aac_mod = 0;
        pac_mod = PAC_mod(i);
        Vhi = Vhi.*(1+pac_mod*s);                       % Modulate high freq activity by modulation envelope.
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


        [XX1,P1] = glmfun(Vlo, Vhi, 'empirical','none',.05);
        rcfc(i,j) = XX1.rcfc; rpac(i,j) = XX1.rpac; rpac_new = XX1.rpac_new;
        p_rcfc(i,j) = P1.rcfc; p_rpac(i,j) = P1.rpac; p_rpac_new = P1.rpac_new;
        [MI1,p_MI1] = modulation_index(Vlo,Vhi,'pvals');
        mi(i,j) = MI1; p_mi(i,j) = p_MI1;
    end
end

strname = ['R_MI_Comparison_Weak_Coupling'];
save(strname)