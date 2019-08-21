%% Code to Generate Quantile Figure
addpath('Chaotic Systems Toolbox')

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

pval = [];
for j = 1:100
    Create_Signals;
    VpinkTest = VpinkTest(2001:end-2000);
    
    amp_LO = 1.2;
    amp_HI = 1;
    V = [ones(1,length(VLOW)/2),amp_LO*ones(1,length(VLOW)/2)].*VLOW + amp_HI*VHI + VpinkTest;

    Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
    Vhi = filtfilt(filtwts_hi,1,V);        

    Vlo_pre = Vlo(1:length(Vlo)/2);
    Vhi_pre = Vhi(1:length(Vhi)/2);
    Vlo_post = Vlo(length(Vlo)/2+1:end);
    Vhi_post = Vhi(length(Vhi)/2+1:end);

    [R,P,I] = glmfun_with_indicator_update(Vlo_pre',Vlo_post',Vhi_pre',Vhi_post','empirical','none',.05);
    pval(j) = P.rpac_condition;
end

strname = ['New_Pval_Condition_Sim_Alow'];
save(strname)