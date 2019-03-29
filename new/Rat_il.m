% Pre
clear;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,4); %four days, 4 electrodes
p_RPAC = zeros(4,4);

RPAC_new = zeros(4,4);
p_RPAC_new = zeros(4,4);

RAAC = zeros(4,4); 
RAAC_new = zeros(4,4);
RCFC = zeros(4,4); 
MI = zeros(4,4);
p_MI = zeros(4,4);

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

%pre
for k = 1:4 %day
        str = ['il_data_pre_day_',num2str(k)];
        load(str)
        for j = 1:4 %electrode
            tic
            j
            il = il_data_pre_day(:,j); 
            il = decimate(il,10);
            il = decimate(il,3);
            
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            Vhi = filtfilt(filtwts_hi,1,il);
            
            [r,p] = glmfun(Vlo',Vhi','empirical','none',.05);
            RPAC_new(k,j) = r.rpac_new; 
            p_RPAC_new(k,j) = p.rpac_new;
            
            RPAC(k,j) = r.rpac;
            p_RPAC(k,j) = p.rpac;
            
            RAAC_new(k,j) = r.raac_new; 
            RAAC(k,j) = r.raac;
            RCFC(k,j) = r.rcfc; 
            
            [mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
            p_MI(k,j) = p_mi;
            MI(k,j) = mi;
            toc
        end
end
save('Rat_Data_Pre_il','RPAC','p_RPAC','RPAC_new','p_RPAC_new','RAAC','RAAC_new','RCFC','MI','p_MI')

% Post
clear;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,4); %four days, 4 electrodes
p_RPAC = zeros(4,4);

RPAC_new = zeros(4,4);
p_RPAC_new = zeros(4,4);

RAAC = zeros(4,4); 
RAAC_new = zeros(4,4);
RCFC = zeros(4,4); 
MI = zeros(4,4);
p_MI = zeros(4,4);

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
        str = ['il_data_post_day_',num2str(k)];
        load(str)
        for j = 1:4 %electrode
            il = il_data_post_day(:,j); 
            il = decimate(il,10);
            il = decimate(il,3);
            
            Vlo = filtfilt(filtwts_lo,1,il);            % Define low freq band activity.
            Vhi = filtfilt(filtwts_hi,1,il);
            
            [r,p] = glmfun(Vlo',Vhi','empirical','none',.05);
            RPAC_new(k,j) = r.rpac_new; 
            p_RPAC_new(k,j) = p.rpac_new;
            
            RPAC(k,j) = r.rpac;
            p_RPAC(k,j) = p.rpac;
            
            RAAC_new(k,j) = r.raac_new; 
            RAAC(k,j) = r.raac;
            RCFC(k,j) = r.rcfc; 
            
            [mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
            p_MI(k,j) = p_mi;
            MI(k,j) = mi;
        end
end
save('Rat_Data_Post_il','RPAC','p_RPAC','RPAC_new','p_RPAC_new','RAAC','RAAC_new','RCFC','MI','p_MI')