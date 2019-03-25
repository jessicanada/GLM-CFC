%% 5th and 95th quantile values from 1000 simulations
load('R_MI_Comparison_Increase_Alow_Res.mat')

% blue is RPAC_new
% red is modulation index
% everything is z-scored

for i = 1:21    
    hold on;
    modulation_level = mod(i);
    
    x1 = RPAC(i,:);
    x2 = RPAC_new(i,:); 
    x3 = MI(i,:); 
    
    sigma1 = std(x1); sigma2 = std(x2); sigma3 = std(x3);
    
    x1 = (x1-mean(RPAC(1,:)))/std(RPAC(1,:));
    x2 = (x2-mean(RPAC_new(1,:)))/std(RPAC_new(1,:));
    x3 = (x3-mean(MI(1,:)))/std(MI(1,:));
    
    mn1 = median(x1);
    lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
    uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
    mn2 = median(x2);
    lq2 = quantile(x2,.95);
    uq2 = quantile(x2,.05);
    mn3 = median(x3);
    lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
    uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
    
%     plot([modulation_level+.01,modulation_level+.01], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
%     plot([modulation_level+.01,modulation_level+.01], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
%     plot(modulation_level+.01, mn1, 'ok', 'MarkerSize', 8)
%     plot(modulation_level+.01, lq1, 'xk')
%     plot(modulation_level+.01, uq1, 'xk')
    plot([modulation_level,modulation_level], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
    plot([modulation_level,modulation_level], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
    plot(modulation_level, mn2, 'ob', 'MarkerSize', 8)
    plot(modulation_level, lq2, 'xb')
    plot(modulation_level, uq2, 'xb')
    plot([modulation_level-.01,modulation_level-.01], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
    plot([modulation_level-.01,modulation_level-.01], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
    plot(modulation_level-.01, mn3, 'or', 'MarkerSize', 8)
    plot(modulation_level-.01, lq3, 'xr')
    plot(modulation_level-.01, uq3, 'xr')
    
    xlim([.8,5.2])
end

%% Confidence Intervals
mod = [1,1.2,1.4,1.6,1.8,2,2.2];
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
Create_Signals;
VpinkTest = VpinkTest(2001:end-2000);
figure(1)
for i = 1:7
        amp_LO = mod(i);
        amp_HI = 5;
        V = amp_LO*VLOW + amp_HI*VHI + VpinkTest;
      
        Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V);        

        nCtlPts = 10;
        [XX] = glmfun(Vlo, Vhi,'none','ci',.05);
        
        
        hold on;
        modulation_level = mod(i);
%         plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(1)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%         plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(2)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%         plot(modulation_level+.02, XX.rpac, 'ok', 'MarkerSize', 8,'LineWidth',2)
%         plot(modulation_level+.02, XX.rpac_ci(1), 'xk','LineWidth',2)
%         plot(modulation_level+.02, XX.rpac_ci(2), 'xk','LineWidth',2)

        plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(1)], 'Color', 'b', 'LineWidth', 2)
        plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(2)], 'Color', 'b', 'LineWidth', 2)
        plot(modulation_level, XX.rpac_new, 'ob', 'MarkerSize', 8,'LineWidth',2)
        plot(modulation_level, XX.rpac_new_ci(1), 'xb','LineWidth',2)
        plot(modulation_level, XX.rpac_new_ci(2), 'xb','LineWidth',2)
        
        
    xlim([.8,2.4])
    title('RPAC (new)')
end

figure(2)
for i = 1:7
        amp_LO = mod(i);
        amp_HI = 5;
        V = amp_LO*VLOW + amp_HI*VHI + VpinkTest;
      
        Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V);        

        nCtlPts = 10;
        [mi] = modulation_index(Vlo, Vhi,'none');
        mi_perm = zeros(1,1000);
        for j = 1:1000
            mi_perm(j) = modulation_index_randperm(Vlo,Vhi);
        end

        mi_ci_low = quantile(mi_perm,.05);
        mi_ci_high = quantile(mi_perm,.95);
        
        hold on;
        modulation_level = mod(i);
        plot([modulation_level+.02,modulation_level+.02], [mi, mi_ci_low], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
        plot([modulation_level+.02,modulation_level+.02], [mi, mi_ci_high], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
        plot(modulation_level+.02, mi, 'ok', 'MarkerSize', 8,'LineWidth',2)
        plot(modulation_level+.02, mi_ci_low, 'xk','LineWidth',2)
        plot(modulation_level+.02, mi_ci_high, 'xk','LineWidth',2)
    
    xlim([.8,2.4])
    title('Modulation Index')
end
% xlabel('Scale Factor of Low-Frequency Amplitude'); ylabel('z-scored R')
% hold off
% set(gca,'FontSize',15)



%% Code to Generate Quantile Figure

mod = [1:.2:5];
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

N = 100;
RPAC = zeros(7,N); RPAC_new = zeros(7,N);RCFC = zeros(7,N);MI = zeros(7,N); RPAC_new_mak = zeros(7,N);

for j = 1:N
    j
    Create_Signals;
    VpinkTest = VpinkTest(2001:end-2000);
    for i = 1:21
        
        amp_LO = mod(i);
        amp_HI = 5;
        V = amp_LO*VLOW + amp_HI*VHI + VpinkTest;

        Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V);        

        nCtlPts = 10;
        [XX] = glmfun(Vlo, Vhi,'none','none',.05);
        RPAC(i,j) = XX.rpac; RPAC_new(i,j) = XX.rpac_new; RCFC(i,j) = XX.rcfc;
        RPAC_new_mak(i,j) = XX.rpac_new_mak;
        MI(i,j) = modulation_index(Vlo,Vhi,'none');
    end
end

strname = ['R_MI_Comparison_Increase_Alow_Res'];
save(strname)
