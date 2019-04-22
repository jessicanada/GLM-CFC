%% Code to Generate Quantile Figure
addpath('Chaotic Systems Toolbox')

mod = [1:.2:2.6];
Create_Signals;

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

N = 1000;
RPAC = zeros(length(mod),N);MI = zeros(length(mod),N);
for j = 1:N
    j
    Create_Signals;
    VpinkTest = VpinkTest(2001:end-2000);
    
    for i = 1:length(mod)
        amp_LO = mod(i);
        amp_HI = 5;
        V = amp_LO*VLOW + amp_HI*VHI + VpinkTest;

        Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V);        

        nCtlPts = 10;
        [XX] = glmfun(Vlo, Vhi,'none','none','none',.05);
        RPAC(i,j) = XX.rpac_new;
        MI(i,j) = modulation_index(Vlo,Vhi,'none');
    end
end

strname = ['Increase_Alow_Results_200_Seconds'];
save(strname,'RPAC','MI')