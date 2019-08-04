% %% PAC AAC
% addpath('Chaotic Systems Toolbox')
% N = 1000;
% 
% %AAC
% RPAC = zeros(1,N);
% RAAC_new = zeros(1,N);
% RPAC_new = zeros(1,N);
% RCFC = zeros(1,N);
% RAAC = zeros(1,N);
% p_PAC = zeros(1,N); p_AAC = zeros(1,N); p_CFC = zeros(1,N);
% p_PAC_new = zeros(1,N); p_AAC_new = zeros(1,N);
% 
% i = str2num(id);
% i
% MOD = [0:.1:1.9];
% mval_AAC = MOD(i);
% 
% for j = 1:N
%     j
%     [XX,P] = simfun(0,mval_AAC,'pink','empirical','none','none',.05);
%     RPAC(j) = XX.rpac;
%     RPAC_new(j) = XX.rpac_new;
%     RCFC(j) = XX.rcfc;
%     RAAC(j) = XX.raac;
%     RAAC_new(j) = XX.raac_new;
%     p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
%     p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
% end
% 
% strname = ['AAC_Simulations_',id];
% save(strname)
% 
% %PAC
% RPAC = zeros(1,N);
% RAAC_new = zeros(1,N);
% RPAC_new = zeros(1,N);
% RCFC = zeros(1,N);
% RAAC = zeros(1,N);
% p_PAC = zeros(1,N); p_AAC = zeros(1,N); p_CFC = zeros(1,N);
% p_PAC_new = zeros(1,N); p_AAC_new = zeros(1,N);
% 
% MOD = [0:.05:.95];
% mval_PAC = MOD(i)
% for j = 1:N
%     [XX,P] = simfun(mval_PAC,0,'pink','empirical','none','none',.05);
%     RPAC(j) = XX.rpac;
%     RPAC_new(j) = XX.rpac_new;
%     RCFC(j) = XX.rcfc;
%     RAAC(j) = XX.raac;
%     RAAC_new(j) = XX.raac_new;
%     p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
%     p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
% end
% 
% strname = ['PAC_Simulations_',id];
% save(strname)
% 
% %CFC
% RPAC = zeros(1,N);
% RAAC_new = zeros(1,N);
% RPAC_new = zeros(1,N);
% RCFC = zeros(1,N);
% RAAC = zeros(1,N);
% p_PAC = zeros(1,N); p_AAC = zeros(1,N); p_CFC = zeros(1,N);
% p_PAC_new = zeros(1,N); p_AAC_new = zeros(1,N);
% 
% 
% for j = 1:N
%     [XX,P] = simfun(mval_PAC,mval_AAC,'pink','empirical','none','none',.05);
%     RPAC(j) = XX.rpac;
%     RPAC_new(j) = XX.rpac_new;
%     RCFC(j) = XX.rcfc;
%     RAAC(j) = XX.raac;
%     RAAC_new(j) = XX.raac_new;
%     p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
%     p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
% end
% 
% strname = ['CFC_Simulations_',id];
% save(strname)
% 
i = str2num(id);
MOD = [0:.025:.5];
mval_PAC = MOD(i);
RAAC_new = zeros(1,N);
RPAC_new = zeros(1,N);
p_PAC_new = zeros(1,N); p_AAC_new = zeros(1,N);
MI = zeros(1,N); p_MI = zeros(1,N);

for j = 1:N
    [XX,P,Vlo,Vhi] = simfun(mval_PAC,0,'pink','empirical','none','none',.05);
    RPAC_new(j) = XX.rpac_new;
    RAAC_new(j) = XX.raac_new;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
    [mi,p] = modulation_index(Vlo,Vhi,'pvals');
    MI(j) = mi; p_MI(j) = p;
end

strname = ['CFC_Simulations_Sensitivity_zero_AAC',id];
save(strname)

% %%
% MOD = [0:.05:.95];
% rpac = zeros(20,1000);
% raac_new = zeros(20,1000);
% rpac_new = zeros(20,1000);
% rcfc = zeros(20,1000);
% raac = zeros(20,1000);
% p_pac = zeros(20,1000); p_aac = zeros(20,1000); p_cfc = zeros(20,1000);
% p_pac_new = zeros(20,1000); p_aac_new = zeros(20,1000);
% for i = 1:20
%     strname = ['PAC_Simulations_',num2str(i)];
%     load(strname)
%     rpac(i,:) = RPAC; rcfc(i,:) = RCFC; raac(i,:) = RAAC; rpac_new(i,:) = RPAC_new; raac_new(i,:) = RAAC_new;
%     p_pac(i,:) = p_PAC; p_aac(i,:) = p_AAC; p_cfc(i,:) = p_CFC; p_pac_new(i,:) = p_PAC_new; p_aac_new(i,:) = p_AAC_new;
% end
% RPAC = rpac; RCFC = rcfc; RAAC = raac; RPAC_new = rpac_new; RAAC_new = raac_new;
% PPAC = p_pac; PAAC = p_aac; PCFC = p_cfc; PPAC_new = p_pac_new; PAAC_new = p_aac_new;
% 
% figure;
% subplot(1,2,1)
% L = 20;
% rpac = RPAC_new(L,:); raac = RAAC_new(L,:); rcfc = RCFC(L,:);
% ppac = PPAC_new(L,:); paac = PAAC_new(L,:); pcfc = PCFC(L,:);
% 
% ind_aac = find(paac<.05); ind_pac = find(ppac<.05); ind_cfc = find(pcfc<.05);
% 
% h1 = histogram(rpac(ind_pac)); hold on; h2 = histogram(rcfc(ind_cfc));
% h3 = histogram(raac(ind_aac));
% legend('PAC','CFC','AAC')
% set(gca,'FontSize',12)
% axis tight
% xlabel('R'); ylabel('Count')
% ylim([0,150])
% 
% subplot(1,2,2)
% 
% hold on;
% for i=1:L
%     x1 = RPAC_new(i,:); ind_pac = find(PPAC_new(i,:)<.05); x1 = x1(ind_pac); %grey
%     x2 = RAAC_new(i,:); ind_aac = find(PAAC_new(i,:)<.05); x2 = x2(ind_aac); %blue
%     x3 = RCFC(i,:); ind_cfc = find(PCFC(i,:)<.05); x3 = x3(ind_cfc); %red
%     
%     min1 = min(x1);  %snval(round(0.25*nensemble));
%     max1 = max(x1);  %snval(round(0.75*nensemble));
%     min2 = min(x2);
%     max2 = max(x2);
%     min3 = min(x3);  %snval(round(0.25*nensemble));
%     max3 = max(x3);  %snval(round(0.75*nensemble));
%     
%     mn1 = median(x1);
%     lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
%     uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
%     mn2 = median(x2);
%     lq2 = quantile(x2,.95);
%     uq2 = quantile(x2,.05);
%     mn3 = median(x3);
%     lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
%     uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
%     
%     
%     modulation_level = MOD(i)*100;
%     
%     plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
%     plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
%     plot(modulation_level+1, mn1, 'ok', 'MarkerSize', 8)
%     plot(modulation_level+1, lq1, 'xk')
%     plot(modulation_level+1, uq1, 'xk')
%     plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
%     plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
%     plot(modulation_level-1, mn2, 'ob', 'MarkerSize', 8)
%     plot(modulation_level-1, lq2, 'xb')
%     plot(modulation_level-1, uq2, 'xb')
%     plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
%     plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
%     plot(modulation_level-1, mn3, 'or', 'MarkerSize', 8)
%     plot(modulation_level-1, lq3, 'xr')
%     plot(modulation_level-1, uq3, 'xr')
%     xlim([0,100])
% end
% set(gca,'FontSize',13)
% grid off
% xlabel('Intensity'); ylabel('R')