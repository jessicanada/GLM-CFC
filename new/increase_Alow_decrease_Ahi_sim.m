diff_pac = zeros(1,100);
diff_mi = zeros(1,100);
for i = 1:100
    dt = 0.002;  Fs = 1/dt;  fNQ = Fs/2;        % Simulated time series parameters.
    N  =20*20/dt+4000;                            % # steps to simulate, making the duration 20s

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
    s = s/max(s);                                       

    aac_mod = [.5*ones(length(Vhi)/2,1);1*ones(length(Vhi)/2,1)]';
    pac_mod = [1*ones(length(Vhi)/2,1);.9*ones(length(Vhi)/2,1)]'; %decrease PAC in post
    Vhi = Vhi.*(1+pac_mod.*s);                       
    Vhi = Vhi.*(1+aac_mod.*AmpLo/max(AmpLo));

    Vpink2 = make_pink_noise(1,N,dt);                   
    noise_level = 0.01;

    Vlo = Vlo.*[1*ones(length(Vhi)/2,1);2*ones(length(Vhi)/2,1)]'; %increase low frequency amplitude in post
    Vhi = Vhi.*[1*ones(length(Vhi)/2,1);.1*ones(length(Vhi)/2,1)]'; %increase high frequency amplitude in post
    V1 = Vlo+Vhi+noise_level*Vpink2;                 	

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

    [XX1] = glmfun(Vlo(1:length(Vhi)/2), Vhi(1:length(Vhi)/2), 'none','none',.05);
    [MI1] = modulation_index(Vlo(1:length(Vhi)/2),Vhi(1:length(Vhi)/2),'none');

    [XX2] = glmfun(Vlo(length(Vhi)/2:end), Vhi(length(Vhi)/2:end), 'none','none',.05);
    [MI2] = modulation_index(Vlo(length(Vhi)/2:end),Vhi(length(Vhi)/2:end),'none');
    
    diff_pac(i) = XX1.rpac-XX2.rpac;
    diff_mi(i) = MI1-MI2;
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(strcat('RPAC: pre = ',num2str(XX1.rpac),' post = ',num2str(XX2.rpac)))
disp(strcat('RCFC: pre = ',num2str(XX1.rcfc),' post = ',num2str(XX2.rcfc)))
disp(strcat('MI: pre = ',num2str(MI1),' post = ',num2str(MI2)))
%%
%histogram(diff_pac/max(diff_pac),10); hold on; histogram(diff_mi/max(diff_mi),10); legend('PAC','MI')
length(find(diff_pac>0 & diff_mi<0))
length(find(diff_mi>0 & diff_pac<0))

%%
[XX,P,I] = glmfun_with_indicator(Vlo(1:5000)',Vlo(5001:end)',Vhi(1:5000)',Vhi(5001:end)','none','none',.05);