K = 100;
pval_nochange_w = zeros(1,K);
pval_nochange_wo = zeros(1,K);
for i = 1:K
    [Vlo_pre,Vhi_pre,Vlo_post,Vhi_post] = simulate_prepost(3,3,0,0);
    [~,~,I] = glmfun_with_indicator_update(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,'none','none',.05);
    pval_nochange_w(i) = I.d_w_pval;
    pval_nochange_wo(i) = I.d_wo_pval;
    i
end

strname = ['New_Pval_Sim'];
save(strname)