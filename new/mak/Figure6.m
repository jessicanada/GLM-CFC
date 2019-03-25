
%%
% Out of 1000 Simulations, how many significant detections do we find?
% look at length of ind variables
% only considering for one value of modulation strength

load('R_MI_Comparison_Amplitude_Dependence_Threshold_75')
ind_pac_dep = find(p_pacdep<.05); %276
ind_pac_new_dep = find(p_pacdep_new<.05); %203
ind_cfc_dep = find(p_cfcdep<.05); %552
ind_pac_ind = find(p_pacind<.05); %183
ind_pac_new_ind = find(p_pacind_new<.05); %91
ind_cfc_ind = find(p_cfcind<.05); %189

ind_mi_dep = find(p_mi_dep<.05); %282
ind_mi_ind = find(p_mi_ind<.05); %213


%% Code to generate above

addpath('Chaotic Systems Toolbox')

N = 10;
p_cfcdep = zeros(1,N); p_pacdep = zeros(1,N);
p_cfcind = zeros(1,N); p_pacind = zeros(1,N);
r_pacdep = zeros(1,N); r_cfcdep = zeros(1,N);
r_pacind = zeros(1,N); r_cfcind = zeros(1,N);
p_mi_dep = zeros(1,N); p_mi_ind = zeros(1,N);

for i = 1:N
    i
    [p,xx,P,XX,mi,p_mi,MI,P_MI] = SimForMark0515(.75,1);
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

strname = ['R_MI_Comparison_Amplitude_Dependence_Threshold_75'];
save(strname)
