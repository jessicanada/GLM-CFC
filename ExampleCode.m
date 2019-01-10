% Note: the following simulations use a GLM-based method to simulate V_low
% and V_high, and the p-values are the analytical p-values from the
% manuscript.

%% Signal with neither PAC nor AAC

[XX,P,Vlo,Vhi,t] = simfun(0,0,'GLM','theoretical','none',.05);
figure(1)
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
legend('V_{low}','V_{high}')
axis off
set(gca,'FontSize',13)
str = strcat({'R_PAC is '},{num2str(XX.rpac)},{' with p-value '},{num2str(P.rpac)}); disp(str);
str = strcat({'R_AAC is '},{num2str(XX.raac)},{' with p-value '},{num2str(P.raac)}); disp(str);
str = strcat({'R_CFC is '},{num2str(XX.rcfc)},{' with p-value '},{num2str(P.rcfc)}); disp(str);
xlim([8,10])

figure(2); clf
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
zlim([.008,.02])

%% Signal with PAC

[XX,P,Vlo,Vhi,t] = simfun(1,0,'GLM','theoretical','none',.05);
figure(1)
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
legend('V_{low}','V_{high}')
axis off
set(gca,'FontSize',13)
str = strcat({'R_PAC is '},{num2str(XX.rpac)},{' with p-value '},{num2str(P.rpac)}); disp(str);
str = strcat({'R_AAC is '},{num2str(XX.raac)},{' with p-value '},{num2str(P.raac)}); disp(str);
str = strcat({'R_CFC is '},{num2str(XX.rcfc)},{' with p-value '},{num2str(P.rcfc)}); disp(str);
xlim([8,10])

figure(2); clf
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
zlim([.008,.02])

%% Signal with AAC

[XX,P,Vlo,Vhi,t] = simfun(0,1,'GLM','theoretical','none',.05);
figure(1)
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
legend('V_{low}','V_{high}')
axis off
set(gca,'FontSize',13)
str = strcat({'R_PAC is '},{num2str(XX.rpac)},{' with p-value '},{num2str(P.rpac)}); disp(str);
str = strcat({'R_AAC is '},{num2str(XX.raac)},{' with p-value '},{num2str(P.raac)}); disp(str);
str = strcat({'R_CFC is '},{num2str(XX.rcfc)},{' with p-value '},{num2str(P.rcfc)}); disp(str);
xlim([8,10])

figure(2); clf
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
zlim([.008,.02])

%% Signal with PAC and AAC

[XX,P,Vlo,Vhi,t] = simfun(1,1,'GLM','theoretical','none',.05);
figure(1)
plot(t,Vlo + .08,t,Vhi,'LineWidth',2); axis off
hold on;
[pks, ipks] = findpeaks(Vlo);
for i = 1:length(ipks)
    ind = ipks(i);
    plot([t(ind),t(ind)],[Vlo(ind)+.07,Vlo(ind)+.09],'r','LineWidth',2)
end
plot([8,8.1],[-.04,-.04],'k','LineWidth',2)
legend('V_{low}','V_{high}')
axis off
set(gca,'FontSize',13)
str = strcat({'R_PAC is '},{num2str(XX.rpac)},{' with p-value '},{num2str(P.rpac)}); disp(str);
str = strcat({'R_AAC is '},{num2str(XX.raac)},{' with p-value '},{num2str(P.raac)}); disp(str);
str = strcat({'R_CFC is '},{num2str(XX.rcfc)},{' with p-value '},{num2str(P.rcfc)}); disp(str);
xlim([8,10])

figure(2); clf
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
zlim([.008,.02])