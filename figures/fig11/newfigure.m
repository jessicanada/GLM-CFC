%%
clear;
load('/Users/jnadalin/Dropbox/Nadalin (1)/Raw_data/Good_phase_locking_N=4/data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,16); %four days, 16 combinations
RAAC = zeros(4,16); 
RCFC = zeros(4,16); 
MI = zeros(4,16);
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
            for i = 1:4 %am trial
                am = amdata{k,2}; %post
                am = am(:,i);
                am = decimate(am,10);
                am = decimate(am,3);
                
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
end
save('results_GLM_post_rsquared','RPAC','RAAC','RCFC','MI')

%%
clear;
load('/Users/jnadalin/Dropbox/Nadalin (1)/Raw_data/Good_phase_locking_N=4/data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,16); %four days, 16 combinations
RAAC = zeros(4,16); 
RCFC = zeros(4,16); 
MI = zeros(4,16);
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
            for i = 1:4 %am trial
                am = amdata{k,1}; %pre
                am = am(:,i);
                am = decimate(am,10);
                am = decimate(am,3);
                
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
end
save('results_GLM_pre_rsquared','RPAC','RAAC','RCFC','MI')

%%
load('data_for_GLM_PAC.mat')
ildata = raw1;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

%post
k = 4;
j = 4; %il trial
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

[r] = glm_rsquared(Vlo',Vhi','none','ci');
rpac_post = r.rpac; rpac_ci_post = r.rpac_ci;
%[mi_pre,p_pre,MI_pre] = modulation_index(Vlo, Vhi,'pvals');

%save('results_GLM_intra_post_ci','RPAC','RAAC','RCFC','MI')
%%


load('data_for_GLM_PAC.mat')
ildata = raw1;
amdata = raw2;
dt = 1e-3;  Fs = 1/dt;  fNQ = Fs/2; 

RPAC = zeros(4,8); %four days, 8 electrodes
RAAC = zeros(4,8); 
RCFC = zeros(4,8); 
MI = zeros(4,8);

k = 4;
j = 4; %il trial
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

[r] = glmfun(Vlo',Vhi','theoretical','ci',.05);
rpac_pre = r.rpac; rpac_ci_pre = r.rpac_ci;

k = 3;
load('results_GLM_pre.mat')
RPAC_pre = RPAC(k,:); RAAC_pre = RAAC(k,:); RCFC_pre = RCFC(k,:); MI_pre = MI(k,:);
load('results_GLM_post_rsquared.mat')
RPAC_post = RPAC(k,:); RAAC_post = RAAC(k,:); RCFC_post = RCFC(k,:); MI_post = MI(k,:);

figure(1);
subplot(2,3,3)
plot([RAAC_pre;RAAC_post],'r','Marker','o','MarkerFaceColor','r'); hold on;
plot([mean(RAAC_pre),mean(RAAC_post)],'r','Marker','o','MarkerFaceColor','r','LineWidth',2,'LineStyle','--')
plot([RPAC_pre;RPAC_post],'b','Marker','o','MarkerFaceColor','b');
plot([mean(RPAC_pre),mean(RPAC_post)],'b','Marker','o','MarkerFaceColor','b','LineWidth',2,'LineStyle','--')
xlim([.95,2.05])
set(gca,'xticklabel',{'pre','','','','','post'})
set(gca,'FontSize',12)
title('R_{PAC}, R_{AAC}')
subplot(2,3,1)
plot([MI_pre;MI_post],'k','Marker','o','MarkerFaceColor','k');
hold on; plot([mean(MI_pre);mean(MI_post)],'k','Marker','o','MarkerFaceColor','k','LineWidth',2,'LineStyle','--');
xlim([.95,2.05])
title('MI')
set(gca,'xticklabel',{'pre','','','','','post'})
set(gca,'FontSize',12)
subplot(2,3,2)
plot([RCFC_pre;RCFC_post],'k','Marker','o','MarkerFaceColor','k');
hold on; plot([mean(RCFC_pre);mean(RCFC_post)],'k','Marker','o','MarkerFaceColor','k','LineWidth',2,'LineStyle','--')
xlim([.95,2.05])
title('R_{CFC}')
set(gca,'xticklabel',{'pre','','','','','post'})
set(gca,'FontSize',12)

load('results_GLM_intra_pre.mat')
MI_pre = MI(4,4);
RPAC_pre = RPAC(4,4);

load('results_GLM_intra_post.mat')
MI_post = MI(4,4);
RPAC_post = RPAC(4,4);


subplot(2,3,4)
ctrs = 1:2;
bar(ctrs,[MI_pre,MI_post])
load('results_MI_ci_post.mat')
% this is significance, not error bars ahhhhhhh
hold on;
MI_ci_post = [quantile(MI,.05),quantile(MI,.95)];
load('results_MI_ci_pre.mat')
MI_ci_pre = [quantile(MI,.05),quantile(MI,.95)];
errorbar(ctrs,[MI_pre,MI_post],[MI_pre-MI_ci_pre(1),MI_post-MI_ci_post(1)],[MI_ci_pre(2)-MI_pre,MI_ci_post(2)-MI_post],'k.')
title('MI')
names={'pre'; 'post'};
set(gca,'xticklabel',names)
set(gca,'FontSize',12)

subplot(2,3,5)
ctrs = 1:2;
bar([RPAC_pre,RPAC_post])
hold on;
errorbar(ctrs,[rpac_pre,rpac_post],[rpac_pre-rpac_ci_pre(1),rpac_post-rpac_ci_post(1)],[rpac_ci_pre(2)-rpac_pre,rpac_ci_post(2)-rpac_post],'k.')
title('R_{PAC}')
names={'pre'; 'post'};
set(gca,'xticklabel',names)
set(gca,'FontSize',12)
