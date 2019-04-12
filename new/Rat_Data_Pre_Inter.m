% Pre
clear;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,4,4); %four days, 4 electrodes
p_RPAC = zeros(4,4,4);

RPAC_new = zeros(4,4,4);
p_RPAC_new = zeros(4,4,4);
p_RAAC_new = zeros(4,4,4);

RAAC = zeros(4,4,4); 
RAAC_new = zeros(4,4,4);
RCFC = zeros(4,4,4); 
MI = zeros(4,4,4);
p_MI = zeros(4,4,4);

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

count = 0

load('data_for_GLM_PAC_decimated.mat')
am_data_pre = raw2(:,1);
il_data_pre = raw1(:,1);
%pre
for k = 1:4 %day
        am_data_pre_day = am_data_pre{k};
        il_data_pre_day = il_data_pre{k};
        
        for j = 1:4 %electrode
            il = il_data_pre_day(:,j);         
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            
            for i = 1:4                
                am = am_data_pre_day(:,i); 
                Vhi = filtfilt(filtwts_hi,1,am);

                [r,p] = glmfun(Vlo',Vhi','empirical','none','none',.05);
                RPAC_new(k,i,j) = r.rpac_new; 
                p_RPAC_new(k,i,j) = p.rpac_new;

                RPAC(k,i,j) = r.rpac;
                p_RPAC(k,i,j) = p.rpac;

                RAAC_new(k,i,j) = r.raac_new; 
                p_RAAC_new(k,i,j) = p.raac_new;
                
                RAAC(k,i,j) = r.raac;
                RCFC(k,i,j) = r.rcfc; 

                [mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
                p_MI(k,i,j) = p_mi;
                MI(k,i,j) = mi;
                count = count+1
            end
        end
end
save('Rat_Data_Pre_Inter','RPAC','p_RPAC','RPAC_new','p_RPAC_new','RAAC','RAAC_new','RCFC','MI','p_MI')
