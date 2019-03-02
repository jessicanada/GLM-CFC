%%
clear;
load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 8 electrodes
RAAC = zeros(4,8); 
RCFC = zeros(4,8); 
MI = zeros(4,8);
%post
for k = 1:4 %day
        count = 1;
        for j = 1:4 %il trial
            il = ildata{k,2}; %post
            il = il(:,j);
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
            
            locutoff = 70;                             % High freq passband = [100, 140] Hz.
            hicutoff = 110;
            filtorder = 10*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vhi = filtfilt(filtwts,1,il);            % Define high freq band activity.
            
            [r] = glm_rsquared(Vlo',Vhi','none','none');
            RPAC(k,count) = r.rpac_new; 
            RAAC(k,count) = r.raac_new; 
            RCFC(k,count) = r.rcfc_new; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1
        end
        for i = 1:4 %am trial
            am = amdata{k,2}; %post
            am = am(:,i);
            am = decimate(am,10);
            am = decimate(am,3);

            locutoff = 5;                               % Low freq passband = [4,7] Hz.
            hicutoff = 8;
            filtorder = 3*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vlo = filtfilt(filtwts,1,am);            % Define low freq band activity.
            
            % Filter into high freq band.
            locutoff = 70;                             % High freq passband = [100, 140] Hz.
            hicutoff = 110;
            filtorder = 10*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vhi = filtfilt(filtwts,1,am);            % Define high freq band activity.

            [r] = glm_rsquared(Vlo',Vhi','none','none');
            RPAC(k,count) = r.rpac_new; 
            RAAC(k,count) = r.raac_new; 
            RCFC(k,count) = r.rcfc_new; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
        
end
save('results_GLM_intra_post','RPAC','RAAC','RCFC','MI')


clear;
load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 8 electrodes
RAAC = zeros(4,8); 
RCFC = zeros(4,8); 
MI = zeros(4,8);
%post
for k = 1:4 %day
        count = 1;
        for j = 1:4 %il trial
            il = ildata{k,1}; %pre
            il = il(:,j);
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
            
            locutoff = 70;                             % High freq passband = [100, 140] Hz.
            hicutoff = 110;
            filtorder = 10*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vhi = filtfilt(filtwts,1,il);            % Define high freq band activity.
            
            [r] = glm_rsquared(Vlo',Vhi','none','none');
            RPAC(k,count) = r.rpac_new; 
            RAAC(k,count) = r.raac_new; 
            RCFC(k,count) = r.rcfc_new; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1
        end
        for i = 1:4 %am trial
            am = amdata{k,1}; %pre
            am = am(:,i);
            am = decimate(am,10);
            am = decimate(am,3);

            locutoff = 5;                               % Low freq passband = [4,7] Hz.
            hicutoff = 8;
            filtorder = 3*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vlo = filtfilt(filtwts,1,am);            % Define low freq band activity.
            
            % Filter into high freq band.
            locutoff = 70;                             % High freq passband = [100, 140] Hz.
            hicutoff = 110;
            filtorder = 10*fix(Fs/locutoff);
            MINFREQ = 0;
            trans          = 0.15;                      % fractional width of transition zones
            f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
            m=[0       0                      1            1            0                      0];
            filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
            Vhi = filtfilt(filtwts,1,am);            % Define high freq band activity.

            [r] = glm_rsquared(Vlo',Vhi','none','none');
            RPAC(k,count) = r.rpac_new; 
            RAAC(k,count) = r.raac_new; 
            RCFC(k,count) = r.rcfc_new; 
            [mi] = modulation_index(Vlo,Vhi,'none');
            MI(k,count) = mi;
            count = count+1;
        end
        
end
save('results_GLM_intra_pre','RPAC','RAAC','RCFC','MI')