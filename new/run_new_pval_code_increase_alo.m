N = 1000;
pval = zeros(1,N);
addpath('Chaotic Systems Toolbox')
for i = 1:N
    i
    [Vlo_pre,Vhi_pre,Vlo_post,Vhi_post] = simulate_prepost_increase_alow(1,1,0,0);
    [XX,P,I] = glmfun_with_indicator_update(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,'empirical','none',.05);
    pval(i) = P.rpac_condition;
end
strname = ['New_Pval_Condition_Sim_Increase_Alow'];
save(strname)