clear;

% Evaluatle R_PAC,
[XX,P,Vlo,Vhi,t] = simfun(0,0,'spiking','empirical','none',.05);
% ... and plot the result.
figure(11); clf; subplot(2,2,1);
surf(XX.ampAXIS,XX.phi0,XX.PAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[35, 106, 185]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.AAC,'EdgeColor','none','FaceAlpha',.8,'FaceColor',[214, 26, 70]/255);
hold on;
surf(XX.ampAXIS,XX.phi0,XX.CFC,'EdgeColor','none','FaceAlpha',.9,'FaceColor',[253, 220, 34]/255);
hold off
xlabel('A_{low}'); ylabel('\Phi_{low}'); zlabel('A_{high}')
set(gca,'FontSize',14)
grid off
title(['Rpacnew = ' num2str(XX.rpac_new,3) ', p = ' num2str(P.rpac_new,3) ')'])

% Evaluate MI,
[mi,p_mi] = modulation_index(Vlo,Vhi,'pvals');
% ... and plot the result.
figure(11); subplot(2,2,2)
phi=angle(hilbert(Vlo));            %Compute phase of low freq signal.
amp=abs(hilbert(Vhi));              %Compute amplitude of high freq signal.
N      = 18;                        %Number of phase bins.
p_bins = linspace(-pi,pi,N+1);       %Define the phase bins.
a_mean = zeros(length(p_bins)-1,1);	%Vector for average amps.
p_mean = zeros(length(p_bins)-1,1);	%Vector for phase bins.
for k=1:N                           %For each phase bin,
    pL = p_bins(k);					%... lower phase limit,
    pR = p_bins(k+1);				%... upper phase limit.
    indices=find(phi>=pL & phi<pR);	%Find phases falling in bin,
    a_mean(k) = mean(amp(indices));	%... compute mean amplitude,
    p_mean(k) = mean([pL, pR]);		%... save center phase.
end
p_normalized = a_mean / sum(a_mean);
bar(p_mean, p_normalized)
xlabel('\Phi_{low}')
ylabel('A_{high}')
title(['MI = ' num2str(mi, 3) ', (p = ' num2str(p_mi,3) ')'])
set(gca,'FontSize',14)

% Plot the signals.
figure(11); subplot(2,2,3:4)
plot(t,Vlo)
hold on
plot(t,Vhi)
hold off
xlabel('Time [s]')
legend({'Vlo', 'Vhi'})