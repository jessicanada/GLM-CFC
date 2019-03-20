
%figure(1);
mod = [1,1.2,1.4,1.6,1.8,2,2.2];

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

RPAC = zeros(1,1000); RPAC_new = zeros(1,1000);RCFC = zeros(1,1000);MI = zeros(1,1000);

i=str2num(id);
for j = 1:1000
        
        Create_Signals;
        VpinkTest = VpinkTest(2001:end-2000);
        amp_LO = mod(i);
        amp_HI = 5;
        V = amp_LO*VLOW + amp_HI*VHI + VpinkTest;
      
        Vlo = filtfilt(filtwts_lo,1,V);            % Define low freq band activity.
        Vhi = filtfilt(filtwts_hi,1,V);        

        nCtlPts = 10;
        [XX] = glmfun(Vlo, Vhi,'none','none',.05);
        RPAC(j) = XX.rpac; RPAC_new(j) = XX.rpac_new; RCFC(j) = XX.rcfc;
        MI(j) = modulation_index(Vlo,Vhi,'none');
        
%     hold on;
%     modulation_level = mod(i);
%     plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(1)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%     plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(2)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%     plot(modulation_level+.02, XX.rpac, 'ok', 'MarkerSize', 8,'LineWidth',2)
%     plot(modulation_level+.02, XX.rpac_ci(1), 'xk','LineWidth',2)
%     plot(modulation_level+.02, XX.rpac_ci(2), 'xk','LineWidth',2)
%     
%     plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(1)], 'Color', 'b', 'LineWidth', 2)
%     plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(2)], 'Color', 'b', 'LineWidth', 2)
%     plot(modulation_level, XX.rpac_new, 'ob', 'MarkerSize', 8,'LineWidth',2)
%     plot(modulation_level, XX.rpac_new_ci(1), 'xb','LineWidth',2)
%     plot(modulation_level, XX.rpac_new_ci(2), 'xb','LineWidth',2)
%     
%     xlim([.8,2.4])
end
% xlabel('Scale Factor of Low-Frequency Amplitude'); ylabel('z-scored R')
% hold off
% set(gca,'FontSize',15)

strname = ['R_MI_Comparison_Increase_Alow_', id];
save(strname)

 %%
% 
% figure(1);
% MOD = 0:.1:1;
% hold on;
% %Blue is rpac_new Black is rpac
% for i=1:length(MOD)
%     [XX] = simfun(MOD(i),0,'pink','none','ci',.05);
%         
%     
%     modulation_level = MOD(i);
%     plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(1)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%     plot([modulation_level+.02,modulation_level+.02], [XX.rpac, XX.rpac_ci(2)], 'Color', [0.8,0.8,0.8], 'LineWidth', 2)
%     plot(modulation_level+.02, XX.rpac, 'ok', 'MarkerSize', 8,'LineWidth',2)
%     plot(modulation_level+.02, XX.rpac_ci(1), 'xk','LineWidth',2)
%     plot(modulation_level+.02, XX.rpac_ci(2), 'xk','LineWidth',2)
%     
%     plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(1)], 'Color', 'b', 'LineWidth', 2)
%     plot([modulation_level,modulation_level], [XX.rpac_new, XX.rpac_new_ci(2)], 'Color', 'b', 'LineWidth', 2)
%     plot(modulation_level, XX.rpac_new, 'ob', 'MarkerSize', 8,'LineWidth',2)
%     plot(modulation_level, XX.rpac_new_ci(1), 'xb','LineWidth',2)
%     plot(modulation_level, XX.rpac_new_ci(2), 'xb','LineWidth',2)
%     
%     xlim([0,1.1])
% end
% xlabel('Scale Factor of Low-Frequency Amplitude'); ylabel('z-scored R')
% hold off
% set(gca,'FontSize',15)