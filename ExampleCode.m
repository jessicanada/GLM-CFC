% Note: the following simulations use a GLM-based method to simulate V_low
% and V_high, and the p-values are the empirical bootstrapped p-values from the
% manuscript

%% Signal with neither PAC nor AAC

[XX] = simfun(0,0,'pink','none','none','none',.05);
d = 5;
XX.ampAXIS = XX.ampAXIS(1:d:end);
XX.phi0 = XX.phi0(1:d:end);
XX.PAC = XX.Phi_low(1:d:end,1:d:end); XX.AAC = XX.A_low(1:d:end,1:d:end); XX.CFC = XX.Phi_low_A_low(1:d:end,1:d:end);

surf(XX.ampAXIS,XX.phi0,XX.PAC,'FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX.ampAXIS,XX.phi0,XX.CFC,'FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('\phi_{low}','A_{low}','\phi_{low},A_{low}')
set(gca,'FontSize',13)
set(gca,'YTick',-pi:pi:pi) 
set(gca,'YTickLabel',{'-\pi','0','\pi'})
set(gca,'Ydir','reverse')
grid off
zlim([.002,.02])

%% Signal with PAC

[XX] = simfun(1,0,'pink','none','none','none',.05);
d = 5;
XX.ampAXIS = XX.ampAXIS(1:d:end);
XX.phi0 = XX.phi0(1:d:end);
XX.PAC = XX.Phi_low(1:d:end,1:d:end); XX.AAC = XX.A_low(1:d:end,1:d:end); XX.CFC = XX.Phi_low_A_low(1:d:end,1:d:end);

surf(XX.ampAXIS,XX.phi0,XX.PAC,'FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX.ampAXIS,XX.phi0,XX.CFC,'FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('\phi_{low}','A_{low}','\phi_{low},A_{low}')
set(gca,'FontSize',13)
set(gca,'YTick',-pi:pi:pi) 
set(gca,'YTickLabel',{'-\pi','0','\pi'})
set(gca,'Ydir','reverse')
grid off
zlim([.002,.01])
%% Signal with AAC

[XX] = simfun(0,1,'pink','none','none','none',.05);
d = 5;
XX.ampAXIS = XX.ampAXIS(1:d:end);
XX.phi0 = XX.phi0(1:d:end);
XX.PAC = XX.Phi_low(1:d:end,1:d:end); XX.AAC = XX.A_low(1:d:end,1:d:end); XX.CFC = XX.Phi_low_A_low(1:d:end,1:d:end);

surf(XX.ampAXIS,XX.phi0,XX.PAC,'FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX.ampAXIS,XX.phi0,XX.CFC,'FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('\phi_{low}','A_{low}','\phi_{low},A_{low}')
set(gca,'FontSize',13)
set(gca,'YTick',-pi:pi:pi) 
set(gca,'YTickLabel',{'-\pi','0','\pi'})
set(gca,'Ydir','reverse')
grid off
zlim([.002,.01])

%% Signal with PAC and AAC

[XX] = simfun(1,1,'pink','none','none','none',.05);
d = 5;
XX.ampAXIS = XX.ampAXIS(1:d:end);
XX.phi0 = XX.phi0(1:d:end);
XX.PAC = XX.Phi_low(1:d:end,1:d:end); XX.AAC = XX.A_low(1:d:end,1:d:end); XX.CFC = XX.Phi_low_A_low(1:d:end,1:d:end);

surf(XX.ampAXIS,XX.phi0,XX.PAC,'FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX.ampAXIS,XX.phi0,XX.CFC,'FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlim([min(XX.ampAXIS),max(XX.ampAXIS)]); ylim([min(XX.phi0),max(XX.phi0)])
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
legend('\phi_{low}','A_{low}','\phi_{low},A_{low}')
set(gca,'FontSize',13)
set(gca,'YTick',-pi:pi:pi) 
set(gca,'YTickLabel',{'-\pi','0','\pi'})
set(gca,'Ydir','reverse')
grid off
zlim([.002,.013])