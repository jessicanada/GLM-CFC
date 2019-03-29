%% Code to analyze single electrode rat data
clear;
load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 8 electrodes
p_RPAC = zeros(4,8);

RPAC_new = zeros(4,8);
p_RPAC_new = zeros(4,8);

RAAC = zeros(4,8); 
RAAC_new = zeros(4,8);
RCFC = zeros(4,8); 
MI = zeros(4,8);
p_MI = zeros(4,8);

locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_lo = firls(filtorder,f,m);             % get FIR filter coefficients

locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_hi = firls(filtorder,f,m);             % get FIR filter coefficients

%post
for k = 1:4 %day
        count = 1;
        for j = 1:4 %il trial
            il = ildata{k,2}; %post
            il = il(:,j);
            il = decimate(il,10);
            il = decimate(il,3);
            
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            Vhi = filtfilt(filtwts_hi,1,il);
            
            [r,p] = glmfun(Vlo',Vhi','empirical','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            p_RPAC_new(k,count) = p.rpac_new;
            
            RPAC(k,count) = r.rpac;
            p_RPAC(k,count) = p.rpac;
            
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            
            [mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
            p_MI(k,count) = p_mi;
            MI(k,count) = mi;
            
            count = count+1;
        end
        for i = 1:4 %am trial
            am = amdata{k,2}; %post
            am = am(:,i);
            am = decimate(am,10);
            am = decimate(am,3);

            % Filter into high freq band.
            Vhi = filtfilt(filtwts_hi,1,am);            % Define high freq band activity.
            Vlo = filtfilt(filtwts_lo,1,am);
            
            [r,p] = glmfun(Vlo',Vhi','empirical','none',.05);
            RPAC_new(k,count) = r.rpac_new; 
            p_RPAC_new(k,count) = p.rpac_new;
            
            RPAC(k,count) = r.rpac;
            p_RPAC(k,count) = p.rpac;
            
            RAAC_new(k,count) = r.raac_new; 
            RAAC(k,count) = r.raac;
            RCFC(k,count) = r.rcfc; 
            
            [mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
            p_MI(k,count) = p_mi;
            MI(k,count) = mi;
            
            count = count+1
        end
end
save('Rat_Data_Post','RPAC','p_RPAC','RPAC_new','p_RPAC_new','RAAC','RAAC_new','RCFC','MI','p_MI')