K = 100;
pval_nochange = zeros(1,K);
for i = 1:K
    [Vlo_pre,Vhi_pre,Vlo_post,Vhi_post] = simulate_prepost(0,0,0,0);
    [~,~,I] = glmfun_with_indicator_update(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,'none','none',.05);
    pval_nochange(i) = I.d_pval;
    i
end

strname = ['New_Pval_Sim'];
save(strname)