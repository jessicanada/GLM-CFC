load('R_MI_Comparison_New_R_Definition_Update.mat')
% sig_MI = # significant detections out of 1000 for MI = 583
% sig_PAC_new = # significant detections out of 1000 for RPAC_new = 958


%% Code to generate above data
addpath('Chaotic Systems Toolbox')

% Example to show where GLM outperforms MI. In this example, we'll generate the 
fixed_val = -1;
%CFC_strength   = (0:.1:1);        % This controls "strength" of PAC.
n_realizations = 10;            % The number of times to re-run each CFC setting.

PAC            = zeros(1, n_realizations);
PAC_new        = zeros(1, n_realizations);
CFC            = zeros(1, n_realizations);
mi             = zeros(1, n_realizations);
p_PAC          = zeros(1, n_realizations);
p_CFC          = zeros(1, n_realizations);
p_mi           = zeros(1, n_realizations);
p_PAC_new      = zeros(1, n_realizations);
p_PAC_new_mak = zeros(1,n_realizations);

%for i_CFC=1:length(CFC_strength)

    for i_run=1:n_realizations
        i_run
        [XX,P,Vlo,Vhi,t] = simfun_mi_vs_r(1.5,fixed_val);
        CFC(i_run) = XX.rcfc;                 % Compute GLM and save it.
        PAC(i_run) = XX.rpac;
        PAC_new(i_run) = XX.rpac_new;
        p_PAC(i_run) = P.rpac;
        p_CFC(i_run) = P.rcfc;
        p_PAC_new(i_run) = P.rpac_new;
        p_PAC_new_mak(i_run) = P.rpac_new_mak;
        [MI,p_MI]  = modulation_index(Vlo, Vhi,'pvals');               % Compute MI and save it.
        mi(i_run) = MI;
        p_mi(i_run) = p_MI;
        
    end



    ind1 = find(p_mi<.05);
    sig_MI = length(ind1);
    ind2 = find(p_CFC<.05);
    sig_CFC = length(ind2);
    ind3 = find(p_PAC<.05);
    sig_PAC = length(ind3);
    ind4 = find(p_PAC_new<.05);
    sig_PAC_new = length(ind4);
    ind5 = find(p_PAC_new_mak<.05);
    sig_PAC_new_mak = length(ind5);
    
%strname = ['R_MI_Comparison_New_R_Definition_Update'];
%save(strname)

