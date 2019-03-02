%%
elec = 4;
load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 
k = 4; %day 4
il = ildata{k,2}; %post
il = il(:,elec); %electrode 4
il = decimate(il,10);
il = decimate(il,3);

locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vlo = filtfilt(filtwts,1,il);            % Define low freq band activity.
Vlo_post = Vlo;
            
locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vhi = filtfilt(filtwts,1,il);            % Define high freq band activity.
Vhi_post = Vhi;

         
[r_post,r_p_post] = glmfun(Vlo',Vhi','empirical','none',.05);
[mi_post,mi_p_post] = modulation_index(Vlo,Vhi,'pvals');


load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 
k = 4; %day 4
il = ildata{k,1}; %pre
il = il(:,elec); %electrode 4
il = decimate(il,10);
il = decimate(il,3);

locutoff = 5;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vlo = filtfilt(filtwts,1,il);            % Define low freq band activity.
Vlo_pre = Vlo;
            
locutoff = 70;                             % High freq passband = [100, 140] Hz.
hicutoff = 110;
filtorder = 10*fix(Fs/locutoff);
MINFREQ = 0;
trans          = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
Vhi = filtfilt(filtwts,1,il);            % Define high freq band activity.
Vhi_pre = Vhi;
           
[r_pre,r_p_pre] = glmfun(Vlo',Vhi','empirical','none',.05);
[mi_pre,mi_p_pre] = modulation_index(Vlo,Vhi,'pvals');

%%
Alo_pre = abs(hilbert(Vlo_pre));
Alo_post = abs(hilbert(Vlo_post));
Ahi_pre = abs(hilbert(Vhi_pre));
Ahi_post = abs(hilbert(Vhi_post));
%%
t = dt*(1:(length(Vlo_pre)+length(Vlo_post)));
t_pre = t(1:length(Vlo_pre));
t_post = t(length(Vlo_pre)+1:end);

plot(t_pre,Vlo_pre,t_pre,Vhi_pre); hold on; plot(t_post,Vlo_post,t_post,Vhi_post)
legend('Vlo pre','Vhi pre','Vlo post','Vhi post')
%%
figure(1)
histogram(Alo_pre,'Normalization','probability'); hold on; histogram(Alo_post,'Normalization','probability'); legend('pre','post')
title('Alo')

figure(2)
histogram(Ahi_pre,'Normalization','probability'); hold on; histogram(Ahi_post,'Normalization','probability'); legend('pre','post')
title('Ahi')
