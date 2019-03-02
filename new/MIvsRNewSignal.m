% Example to show where GLM outperforms MI. In this example, we'll generate the 
fixed_val = -1;
CFC_strength   = (0:.1:1);        % This controls "strength" of PAC.
n_realizations = 100;            % The number of times to re-run each CFC setting.

PAC            = zeros(length(CFC_strength), n_realizations);
PAC_new        = zeros(length(CFC_strength), n_realizations);
CFC            = zeros(length(CFC_strength), n_realizations);
mi             = zeros(length(CFC_strength), n_realizations);
p_PAC          = zeros(length(CFC_strength), n_realizations);
p_CFC          = zeros(length(CFC_strength), n_realizations);
p_mi           = zeros(length(CFC_strength), n_realizations);
p_PAC_new      = zeros(length(CFC_strength), n_realizations);

for i_CFC=1:length(CFC_strength)
i_CFC
    for i_run=1:n_realizations

        [XX,P,Vlo,Vhi,t] = simfun_mi_vs_r(CFC_strength(i_CFC),fixed_val);
        CFC(i_CFC,i_run) = XX.rcfc;                 % Compute GLM and save it.
        PAC(i_CFC,i_run) = XX.rpac;
        PAC_new(i_CFC,i_run) = XX.rpac_new;
        p_PAC(i_CFC,i_run) = P.rpac;
        p_CFC(i_CFC,i_run) = P.rcfc;
        p_PAC_new(i_CFC,i_run) = P.rpac_new;
        [MI,p_MI]  = modulation_index(Vlo, Vhi,'pvals');               % Compute MI and save it.
        mi(i_CFC,i_run) = MI;
        p_mi(i_CFC,i_run) = p_MI;
        
    end
    
end

strname = ['R_MI_Comparison_New_R_Definition'];
save(strname)

% %%
% %load('R_MI_Comparison_Empirical')
% sig_MI = zeros(1,11); sig_CFC = zeros(1,11); sig_PAC = zeros(1,11);
% for i = 1:11
%     ind1 = find(p_mi(i,:)<.05);
%     sig_MI(i) = length(ind1);
%     ind2 = find(p_CFC(i,:)<.05);
%     sig_CFC(i) = length(ind2);
%     ind3 = find(p_PAC(i,:)<.05);
%     sig_PAC(i) = length(ind3);
%     ind4 = find(p_PAC_new(i,:)<.05);
%     sig_PAC_new(i) = length(ind4);
% end
% %%
% plot(CFC_strength,sig_MI,CFC_strength,sig_CFC,CFC_strength,sig_PAC,'LineWidth',2); legend('MI','CFC','PAC')
% set(gca,'FontSize',14)
% xlabel('Intensity'); ylabel('Significant Detections')
% %%
% CFC_strength   = (0:.1:1);
% x_axis = CFC_strength;
% figure;
% subplot(2,1,1)
% plot(x_axis,mi_mean,'b','LineWidth',2)
% hold on
% plot(x_axis,mi_lo, 'b:', x_axis,mi_hi, 'b:','LineWidth',2)
% hold off
% ylabel('MI')
% xlim([0,1])
% set(gca,'FontSize',15)
% 
% subplot(2,1,2)
% plot(x_axis,CFC_mean,'r','LineWidth',2)
% hold on
% plot(x_axis,CFC_lo, 'r:', x_axis,CFC_hi, 'r:','LineWidth',2)
% xlim([0,1])
% hold off
% 
% xlabel('Intensity'); ylabel('R_{CFC}')
% set(gca,'FontSize',15)

