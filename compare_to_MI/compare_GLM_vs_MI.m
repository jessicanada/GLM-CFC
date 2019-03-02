% Example to show where GLM outperforms MI. In this example, we'll generate the 
fixed_val = -1;
CFC_strength   = (0:.1:1);        % This controls "strength" of PAC.
n_realizations = 1000;            % The number of times to re-run each CFC setting.

PAC            = zeros(length(CFC_strength), n_realizations);
mi             = zeros(length(CFC_strength), n_realizations);
p_PAC          = zeros(length(CFC_strength), n_realizations);
p_CFC          = zeros(length(CFC_strength), n_realizations);
p_mi           = zeros(length(CFC_strength), n_realizations);

for i_CFC=1:length(CFC_strength)

    for i_run=1:n_realizations

        [XX,P,Vlo,Vhi,t] = new_simfun(CFC_strength(i_CFC),fixed_val);
        PAC(i_CFC,i_run) = XX.rcfc;                 % Compute GLM and save it.
        p_PAC(i_CFC,i_run) = P.rpac;
        p_CFC(i_CFC,i_run) = P.rcfc;
        [MI,p_MI]  = modulation_index(Vlo, Vhi,'pvals');               % Compute MI and save it.
        mi(i_CFC,i_run) = MI;
        p_mi(i_CFC,i_run) = p_MI;
        
    end
    
end

strname = ['CFC_MI_Comparison_Theoretical'];
save(strname)

%%
load('CFC_MI_Comparison_Theoretical.mat')
%load('CFC_MI_Comparison_Empirical.mat')

ind1 = find(p_mi<.05);
mi_mean = zeros(1,11); mi_lo = zeros(1,11); mi_hi = zeros(1,11);
CFC_mean = zeros(1,11); CFC_lo = zeros(1,11); CFC_hi = zeros(1,11);
sig_MI = zeros(1,11); sig_CFC = zeros(1,11);
for i = 1:11
    MI = mi(i,:);
    ind1 = find(p_mi(i,:)<.05);
    sig_MI(i) = length(ind1);
    MI = MI(ind1);
    mi_mean(i) = mean(MI); mi_lo(i) = quantile(MI,0.025); mi_hi(i) = quantile(MI,0.975);
    CFC = PAC(i,:);
    ind2 = find(p_CFC(i,:)<.05);
    CFC = CFC(ind2);
    sig_CFC(i) = length(ind2);
    CFC_mean(i) = mean(CFC); CFC_lo(i) = quantile(CFC,0.025); CFC_hi(i) = quantile(CFC,0.975);
end
%%
plot(CFC_strength,sig_MI,CFC_strength,sig_CFC,'LineWidth',2); legend('MI','CFC')
set(gca,'FontSize',14)
xlabel('Intensity'); ylabel('Significant Detections')
%%
CFC_strength   = (0:.1:1);
x_axis = CFC_strength;
figure;
subplot(2,1,1)
plot(x_axis,mi_mean,'b','LineWidth',2)
hold on
plot(x_axis,mi_lo, 'b:', x_axis,mi_hi, 'b:','LineWidth',2)
hold off
ylabel('MI')
xlim([0,1])
set(gca,'FontSize',15)

subplot(2,1,2)
plot(x_axis,CFC_mean,'r','LineWidth',2)
hold on
plot(x_axis,CFC_lo, 'r:', x_axis,CFC_hi, 'r:','LineWidth',2)
xlim([0,1])
hold off

xlabel('Intensity'); ylabel('R_{CFC}')
set(gca,'FontSize',15)

