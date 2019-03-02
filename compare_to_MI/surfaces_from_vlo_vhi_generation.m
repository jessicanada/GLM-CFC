%A simulation to try and figure out why Mark's code produces a smooth
%result and mine does not

%create the signals, CFC_strength somewhere between 0 and 20
CFC_strength = 20;
[Vlo, Vhi, t, Alo] = generate_interesting_vlo_and_vhi(CFC_strength);

[R,CHI,XX,bC] = glmfun_test(Vlo, Vhi', 10,'noplot');
%glmfun_test evaluates only at AmpLo = exp(bC)

%%
figure(1)
surf(1:600,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(1:600,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(1:600,XX.phi0,XX.PACAAC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
%plot3([exp(bC),exp(bC)],[-3.14,3.14],[0,0],'LineWidth',2) %this line shows the value glmfun_test evaluates at
legend('PAC','AAC','CFC')
xlabel('A_{lo} = Constant'); ylabel('\Phi_{lo}'); zlabel('A_{hi}')
set(gca,'FontSize',13)
title('New')


%%
[R,CHI,XX] = glmfun(Vlo, Vhi', 10,'noplot');
%glmfun evaluates at all values of AmpLo
figure(2)
%%

%surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);

%surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
surf(XX.ampAXIS,XX.phi0,XX.null,'EdgeColor','none','FaceAlpha',.8,'FaceColor','k');
hold on;
surf(XX.ampAXIS,XX.phi0,XX.PACAAC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[35, 106, 185]/255);
xlabel('A_{lo}'); ylabel('\Phi_{lo}'); zlabel('A_{hi}')
set(gca,'FontSize',13)
[m,I] = max(abs(1-XX.PACAAC./XX.null)); [~,j] = max(m);
plot3([XX.ampAXIS(j),XX.ampAXIS(j)],[XX.phi0(I(j)),XX.phi0(I(j))],[XX.null(I(j),(j)),XX.PACAAC(I(j),(j))],'r','LineWidth',2);
%plot3([exp(bC),exp(bC)],[-3.14,3.14],[0,0],'LineWidth',2) %this line shows the value glmfun_test evaluates at
%legend('PAC','AAC','CFC')
%title('Old')

%%
[R,CHI,XX] = glmfun(Vlo, VHI', 10,'noplot');
%glmfun evaluates at all values of AmpLo
figure(2)
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.PACAAC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlabel('A_{lo}'); ylabel('\Phi_{lo}'); zlabel('A_{hi}')
set(gca,'FontSize',13)
%plot3([exp(bC),exp(bC)],[-3.14,3.14],[0,0],'LineWidth',2) %this line shows the value glmfun_test evaluates at
legend('PAC','AAC','CFC')
title('Old')

figure(3)
[R,CHI,XX] = glmfun(Vlo, VHI2', 10,'noplot');
%glmfun evaluates at all values of AmpLo
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.PACAAC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
xlabel('A_{lo}'); ylabel('\Phi_{lo}'); zlabel('A_{hi}')
set(gca,'FontSize',13)
%plot3([exp(bC),exp(bC)],[-3.14,3.14],[0,0],'LineWidth',2) %this line shows the value glmfun_test evaluates at
legend('PAC','AAC','CFC')
title('Using a different s')

