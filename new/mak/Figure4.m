%% PAC 0 AAC 0
MOD = [0:.05:.95];
[XX,P,Vlo,Vhi,t] = simfun(0,0,'pink','none',.05);

figure(1);
set(gcf, 'Units','centimeters','Position',[0,0,50,50])
subplot(4,4,1)
plot(t,Vlo+.08,t,Vhi,'LineWidth',2);
xlim([8,10])
hold on;
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
axis off
legend('Vlo','Vhi')
set(gca,'FontSize',13)

subplot(4,4,2)
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('PAC','AAC','CFC')
set(gca,'FontSize',13)
grid off
zlim([.002,.02])

subplot(4,4,3)
load('PAC_Simulations_Res.mat')
rpac = RPAC_new(1,:); raac = RAAC_new(1,:); rcfc = RCFC(1,:);
ppac = PPAC_new(1,:); paac = PAAC_new(1,:); pcfc = PCFC(1,:);

ind_aac = find(paac<.05); ind_pac = find(ppac<.05); ind_cfc = find(pcfc<.05);
h1 = histogram(rpac(ind_pac)); hold on; h2 = histogram(rcfc(ind_cfc));
h3 = histogram(raac(ind_aac));
ylim([0,150]) 
set(gca,'FontSize',12)
xlabel('R'); ylabel('Count')

% PAC 1 AAC 0
[XX,P,Vlo,Vhi,t] = simfun(1,0,'pink','none',.05);

subplot(4,4,5)
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
xlim([8,10])
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
axis off

subplot(4,4,6)
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
set(gca,'FontSize',13)
grid off
zlim([.004,.007])

subplot(4,4,7)
load('PAC_Simulations_Res.mat')
rpac = RPAC_new(20,:); raac = RAAC_new(20,:); rcfc = RCFC(20,:);
ppac = PPAC_new(20,:); paac = PAAC_new(20,:); pcfc = PCFC(20,:);

ind_aac = find(paac<.05); ind_pac = find(ppac<.05); ind_cfc = find(pcfc<.05);

h1 = histogram(rpac(ind_pac)); hold on; h2 = histogram(rcfc(ind_cfc));
h3 = histogram(raac(ind_aac));
legend('PAC','CFC','AAC')
set(gca,'FontSize',12)
axis tight
xlabel('R'); ylabel('Count')
ylim([0,150])

subplot(4,4,8)
hold on;
for i=1:20
    x1 = RPAC_new(i,:); ind_pac = find(PPAC_new(i,:)<.05); 
    x1 = x1(ind_pac);
    x2 = RAAC_new(i,:); ind_aac = find(PAAC_new(i,:)<.05); 
    x2 = x2(ind_aac);
    x3 = RCFC(i,:); ind_cfc = find(PCFC(i,:)<.05); 
    x3 = x3(ind_cfc);
    
    min1 = min(x1);  %snval(round(0.25*nensemble));
    max1 = max(x1);  %snval(round(0.75*nensemble));
    min2 = min(x2);
    max2 = max(x2);
    min3 = min(x3);  %snval(round(0.25*nensemble));
    max3 = max(x3);  %snval(round(0.75*nensemble));
    
    mn1 = median(x1);
    lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
    uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
    mn2 = median(x2);
    lq2 = quantile(x2,.95);
    uq2 = quantile(x2,.05);
    mn3 = median(x3);
    lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
    uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
    
    
    modulation_level = MOD(i)*100;
    if length(x1)<=50
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'o', 'MarkerSize', 8,'Color',[0.8,0.8,0.8])
        plot(modulation_level+1, lq1, 'x','Color',[0.8,0.8,0.8])
        plot(modulation_level+1, uq1, 'x','Color',[0.8,0.8,0.8])
    else
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', 'k', 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', 'k', 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'ok', 'MarkerSize', 8)
        plot(modulation_level+1, lq1, 'xk')
        plot(modulation_level+1, uq1, 'xk')
    end
    if length(x2)<=50
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'o', 'MarkerSize', 8, 'Color', [140,184,255]/255)
        plot(modulation_level-1, lq2, 'x', 'Color', [140,184,255]/255)
        plot(modulation_level-1, uq2, 'x', 'Color', [140,184,255]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'ob', 'MarkerSize', 8)
        plot(modulation_level-1, lq2, 'xb')
        plot(modulation_level-1, uq2, 'xb')
    end
    if length(x3)<=50
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'o', 'MarkerSize', 8, 'Color', [255,133,112]/255)
        plot(modulation_level-1, lq3, 'x', 'Color', [255,133,112]/255)
        plot(modulation_level-1, uq3, 'x', 'Color', [255,133,112]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'or', 'MarkerSize', 8)
        plot(modulation_level-1, lq3, 'xr')
        plot(modulation_level-1, uq3, 'xr')
    end
    xlim([0,100])
end
set(gca,'FontSize',13)
grid off
xlabel('Intensity'); ylabel('R')

% PAC 0 AAC 1
subplot(4,4,9)
[XX,P,Vlo,Vhi,t] = simfun(0,1,'pink','theoretical',.05);
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
axis off
xlim([8,10])

subplot(4,4,10)
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);

xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
set(gca,'FontSize',13)
grid off
zlim([.004,.01])

load('AAC_Simulations_Res.mat')
rpac = RPAC_new(20,:); raac = RAAC_new(20,:); rcfc = RCFC(20,:);
ppac = PPAC_new(20,:); paac = PAAC_new(20,:); pcfc = PCFC(20,:);

ind_aac = find(paac<.05); ind_pac = find(ppac<.05); ind_cfc = find(pcfc<.05);
subplot(4,4,11)
h1 = histogram(rpac(ind_pac)); hold on; h2 = histogram(rcfc(ind_cfc));
h3 = histogram(raac(ind_aac));
set(gca,'FontSize',12)
axis tight
xlabel('R'); ylabel('Count')
ylim([0,150])


subplot(4,4,12)
hold on;
for i=1:20
    x1 = RPAC_new(i,:); ind_pac = find(PPAC_new(i,:)<.05); 
    x1 = x1(ind_pac);
    x2 = RAAC_new(i,:); ind_aac = find(PAAC_new(i,:)<.05); 
    x2 = x2(ind_aac);
    x3 = RCFC(i,:); ind_cfc = find(PCFC(i,:)<.05); 
    x3 = x3(ind_cfc);
    
    min1 = min(x1);  %snval(round(0.25*nensemble));
    max1 = max(x1);  %snval(round(0.75*nensemble));
    min2 = min(x2);
    max2 = max(x2);
    min3 = min(x3);  %snval(round(0.25*nensemble));
    max3 = max(x3);  %snval(round(0.75*nensemble));
    
    mn1 = median(x1);
    lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
    uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
    mn2 = median(x2);
    lq2 = quantile(x2,.95);
    uq2 = quantile(x2,.05);
    mn3 = median(x3);
    lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
    uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
    
    
    modulation_level = MOD(i)*100;
    if length(x1)<=50
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'o', 'MarkerSize', 8,'Color',[0.8,0.8,0.8])
        plot(modulation_level+1, lq1, 'x','Color',[0.8,0.8,0.8])
        plot(modulation_level+1, uq1, 'x','Color',[0.8,0.8,0.8])
    else
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', 'k', 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', 'k', 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'ok', 'MarkerSize', 8)
        plot(modulation_level+1, lq1, 'xk')
        plot(modulation_level+1, uq1, 'xk')
    end
    if length(x2)<=50
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'o', 'MarkerSize', 8, 'Color', [140,184,255]/255)
        plot(modulation_level-1, lq2, 'x', 'Color', [140,184,255]/255)
        plot(modulation_level-1, uq2, 'x', 'Color', [140,184,255]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'ob', 'MarkerSize', 8)
        plot(modulation_level-1, lq2, 'xb')
        plot(modulation_level-1, uq2, 'xb')
    end
    if length(x3)<=50
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'o', 'MarkerSize', 8, 'Color', [255,133,112]/255)
        plot(modulation_level-1, lq3, 'x', 'Color', [255,133,112]/255)
        plot(modulation_level-1, uq3, 'x', 'Color', [255,133,112]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'or', 'MarkerSize', 8)
        plot(modulation_level-1, lq3, 'xr')
        plot(modulation_level-1, uq3, 'xr')
    end
    xlim([0,100])
end
set(gca,'FontSize',13)
grid off
xlabel('Intensity'); ylabel('R')

% PAC = 1 AAC = 1 (PAC and AAC)
subplot(4,4,13)
[XX,P,Vlo,Vhi,t] = simfun(1,1,'pink','theoretical',.05);
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
axis off
xlim([8,10])

subplot(4,4,14)
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);

xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
set(gca,'FontSize',13)
grid off
zlim([.003,.015])

subplot(4,4,15)
load('CFC_Simulations_Res.mat')
rpac = RPAC_new(20,:); raac = RAAC_new(20,:); rcfc = RCFC(20,:);
ppac = PPAC_new(20,:); paac = PAAC_new(20,:); pcfc = PCFC(20,:);

ind_aac = find(paac<.05); ind_pac = find(ppac<.05); ind_cfc = find(pcfc<.05);

h1 = histogram(rpac(ind_pac)); hold on; h2 = histogram(rcfc(ind_cfc));
h3 = histogram(raac(ind_aac));

set(gca,'FontSize',12)
axis tight
xlabel('R'); ylabel('Count')
ylim([0,150])


subplot(4,4,16)
hold on;
for i=1:20
    x1 = RPAC_new(i,:); ind_pac = find(PPAC_new(i,:)<.05); 
    x1 = x1(ind_pac);
    x2 = RAAC_new(i,:); ind_aac = find(PAAC_new(i,:)<.05); 
    x2 = x2(ind_aac);
    x3 = RCFC(i,:); ind_cfc = find(PCFC(i,:)<.05); 
    x3 = x3(ind_cfc);
    
    min1 = min(x1);  %snval(round(0.25*nensemble));
    max1 = max(x1);  %snval(round(0.75*nensemble));
    min2 = min(x2);
    max2 = max(x2);
    min3 = min(x3);  %snval(round(0.25*nensemble));
    max3 = max(x3);  %snval(round(0.75*nensemble));
    
    mn1 = median(x1);
    lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
    uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
    mn2 = median(x2);
    lq2 = quantile(x2,.95);
    uq2 = quantile(x2,.05);
    mn3 = median(x3);
    lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
    uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
    
    
    modulation_level = MOD(i)*100;
    if length(x1)<=50
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'o', 'MarkerSize', 8,'Color',[0.8,0.8,0.8])
        plot(modulation_level+1, lq1, 'x','Color',[0.8,0.8,0.8])
        plot(modulation_level+1, uq1, 'x','Color',[0.8,0.8,0.8])
    else
        plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', 'k', 'LineWidth', 1)
        plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', 'k', 'LineWidth', 1)
        plot(modulation_level+1, mn1, 'ok', 'MarkerSize', 8)
        plot(modulation_level+1, lq1, 'xk')
        plot(modulation_level+1, uq1, 'xk')
    end
    if length(x2)<=50
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', [140,184,255]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'o', 'MarkerSize', 8, 'Color', [140,184,255]/255)
        plot(modulation_level-1, lq2, 'x', 'Color', [140,184,255]/255)
        plot(modulation_level-1, uq2, 'x', 'Color', [140,184,255]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
        plot(modulation_level-1, mn2, 'ob', 'MarkerSize', 8)
        plot(modulation_level-1, lq2, 'xb')
        plot(modulation_level-1, uq2, 'xb')
    end
    if length(x3)<=50
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', [255,133,112]/255, 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'o', 'MarkerSize', 8, 'Color', [255,133,112]/255)
        plot(modulation_level-1, lq3, 'x', 'Color', [255,133,112]/255)
        plot(modulation_level-1, uq3, 'x', 'Color', [255,133,112]/255)
    else
        plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
        plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
        plot(modulation_level-1, mn3, 'or', 'MarkerSize', 8)
        plot(modulation_level-1, lq3, 'xr')
        plot(modulation_level-1, uq3, 'xr')
    end
    xlim([0,100])
end
ylim([.13,1.1])
set(gca,'FontSize',13)
grid off
xlabel('Intensity'); ylabel('R')

%%
load('PAC_Simulations_Res.mat')
hold on;
for i=1:20
    x1 = RPAC_new(i,:); ind_pac = find(PPAC_new(i,:)<.05); 
    if length(ind_pac)<=5
        ind_pac = [];
    end
    x1 = x1(ind_pac);
    x2 = RAAC_new(i,:); ind_aac = find(PAAC_new(i,:)<.05); 
    if length(ind_aac)<=5
        ind_aac = [];
    end
    x2 = x2(ind_aac);
    x3 = RCFC(i,:); ind_cfc = find(PCFC(i,:)<.05); 
    if length(ind_cfc)<=5
        ind_cfc = [];
    end
    x3 = x3(ind_cfc);
    
    min1 = min(x1);  %snval(round(0.25*nensemble));
    max1 = max(x1);  %snval(round(0.75*nensemble));
    min2 = min(x2);
    max2 = max(x2);
    min3 = min(x3);  %snval(round(0.25*nensemble));
    max3 = max(x3);  %snval(round(0.75*nensemble));
    
    mn1 = median(x1);
    lq1 = quantile(x1,.95);  %snval(round(0.25*nensemble));
    uq1 = quantile(x1,.05);  %snval(round(0.75*nensemble));
    mn2 = median(x2);
    lq2 = quantile(x2,.95);
    uq2 = quantile(x2,.05);
    mn3 = median(x3);
    lq3 = quantile(x3,.95);  %snval(round(0.25*nensemble));
    uq3 = quantile(x3,.05);  %snval(round(0.75*nensemble));
    
    
    modulation_level = MOD(i)*100;
    
    plot([modulation_level+1,modulation_level+1], [mn1, lq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
    plot([modulation_level+1,modulation_level+1], [mn1, uq1], 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
    plot(modulation_level+1, mn1, 'ok', 'MarkerSize', 8)
    plot(modulation_level+1, lq1, 'xk')
    plot(modulation_level+1, uq1, 'xk')
    plot([modulation_level-1,modulation_level-1], [mn2, lq2], 'Color', 'b', 'LineWidth', 1)
    plot([modulation_level-1,modulation_level-1], [mn2, uq2], 'Color', 'b', 'LineWidth', 1)
    plot(modulation_level-1, mn2, 'ob', 'MarkerSize', 8)
    plot(modulation_level-1, lq2, 'xb')
    plot(modulation_level-1, uq2, 'xb')
    plot([modulation_level-1,modulation_level-1], [mn3, lq3], 'Color', 'r', 'LineWidth', 1)
    plot([modulation_level-1,modulation_level-1], [mn3, uq3], 'Color', 'r', 'LineWidth', 1)
    plot(modulation_level-1, mn3, 'or', 'MarkerSize', 8)
    plot(modulation_level-1, lq3, 'xr')
    plot(modulation_level-1, uq3, 'xr')
    xlim([0,100])
end

%% Code to generate above figures

addpath('Chaotic Systems Toolbox')

%AAC
MOD = [0:.05:.95];
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);

i = str2num(id);
i
mval = MOD(i);

for j = 1:1000
    j
    [XX,P] = simfun(0,mval,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['AAC_Simulations_',id];
save(strname)

%PAC
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);


for j = 1:1000
    [XX,P] = simfun(mval,0,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['PAC_Simulations_',id];
save(strname)

%CFC
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);


for j = 1:1000
    [XX,P] = simfun(mval,mval,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['CFC_Simulations_',id];
save(strname)