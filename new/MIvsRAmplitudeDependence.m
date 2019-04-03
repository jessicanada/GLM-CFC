addpath('Chaotic Systems Toolbox')

N = 1000;
p_cfcdep = zeros(1,N); p_pacdep = zeros(1,N);
p_cfcind = zeros(1,N); p_pacind = zeros(1,N);
r_pacdep = zeros(1,N); r_cfcdep = zeros(1,N);
r_pacind = zeros(1,N); r_cfcind = zeros(1,N);
mi_dep = zeros(1,N); mi_ind = zeros(1,N);
p_mi_dep = zeros(1,N); p_mi_ind = zeros(1,N);
r_pacdep_new = zeros(1,N); r_pacind_new = zeros(1,N);
p_pacdep_new = zeros(1,N); p_pacind_new = zeros(1,N);

for i = 1:N
    i
    [p,xx,P,XX,mi,p_mi,MI,P_MI] = ind_dep_sim(.95,5);
    %simulations where PAC events dependent on low-freq amplitude
    r_pacdep(i) = xx.rpac; %R values
    r_pacdep_new(i) = xx.rpac_new;
    r_cfcdep(i) = xx.rcfc;
    p_pacdep(i) = p.rpac;
    p_pacdep_new(i) = p.rpac_new;
    p_cfcdep(i) = p.rcfc;
    mi_dep(i) = mi;
    p_mi_dep(i) = p_mi;
    %simulations where PAC events independent of low-freq amplitude
    r_pacind(i) = XX.rpac; %R values
    r_pacind_new(i) = XX.rpac_new;
    r_cfcind(i) = XX.rcfc;
    p_pacind(i) = P.rpac;
    p_pacind_new(i) = P.rpac_new;
    p_cfcind(i) = P.rcfc;
    mi_ind(i) = MI;
    p_mi_ind(i) = P_MI;
end

strname = ['R_MI_Comparison_Amplitude_Dependence_Threshold_95_Mod_5'];
save(strname)

 %%
% %load('R_MI_Comparison_Amplitude_Dependence_Threshold_75')
% load('R_MI_Comparison_Amplitude_Dependence_Threshold_90_Mod_2.mat')
% %p_cfcdep = p_cfcdep(1:100); p_cfcind = p_cfcind(1:100);
% ind_pac_dep = find(p_pacdep<.05); %83
% ind_pac_new_dep = find(p_pacdep_new<.05); %27
% ind_cfc_dep = find(p_cfcdep<.05); %48
% ind_pac_ind = find(p_pacind<.05); %102
% ind_pac_new_ind = find(p_pacind_new<.05); %23
% ind_cfc_ind = find(p_cfcind<.05); %32
% 
% ind_mi_dep = find(p_mi_dep<.05); %85
% ind_mi_ind = find(p_mi_ind<.05); %100

% %%
% figure(1)
% %rcfc_dep = r_cfcdep(1:100); r_cfcind = r_cfcind(1:100);
% histogram(r_cfcdep(ind_cfc_dep)); hold on; histogram(r_cfcind(ind_cfc_ind)); legend('dependent','independent'); title('CFC')
% 
% figure(2)
%histogram(mi_dep(ind_mi_dep)); hold on; histogram(mi_ind(ind_mi_ind)); legend('dependent','independent'); title('MI')