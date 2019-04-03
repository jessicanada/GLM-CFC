load('MG49S45_vars.mat')
Ahi = abs(hilbert(Vhi));
Alo = abs(hilbert(Vlo));

%% Visualize fits

figure;
subplot(2,1,1)
histfit(Ahi,2000,'gamma')
title('Ahi')
subplot(2,1,2)
histfit(Alo,2000,'gamma')
title('Alo')

%% compare eCDFs
figure
subplot(2,1,1)
[f,x_values] = ecdf(Ahi);
J = plot(x_values,f,'LineWidth',2);
hold on;
phat_hi = gamfit(Ahi);
K = plot(x_values,gamcdf(x_values,phat_hi(1),phat_hi(2)),'r--','LineWidth',2);
title('Ahi')

subplot(2,1,2)
[f,x_values] = ecdf(Alo);
J = plot(x_values,f,'LineWidth',2);
hold on;
phat_lo = gamfit(Alo);
K = plot(x_values,gamcdf(x_values,phat_lo(1),phat_lo(2)),'r--','LineWidth',2);
title('Alo')

%% K-S test of fit

%reject the null hypothesis that this data does not come from the specified
%distribution

[f,x_values] = ecdf(Ahi);
kstest(Ahi,'CDF',[x_values,gamcdf(x_values,phat_hi(1),phat_hi(2))])


[f,x_values] = ecdf(Alo);
kstest(Alo,'CDF',[x_values,gamcdf(x_values,phat_lo(1),phat_lo(2))])
