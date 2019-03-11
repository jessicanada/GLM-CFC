%% PAC AAC
addpath('Chaotic Systems Toolbox')

%AAC
MOD = [0:.05:.95];
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);

i = str2num(id);
i
mval = MOD(i);

for j = 1:1000
    j
    [XX,P] = simfun(0,mval,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['AAC_Simulations_',id];
save(strname)

%PAC
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);


for j = 1:1000
    [XX] = simfun(mval,0,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['PAC_Simulations_',id];
save(strname)

%CFC
RPAC = zeros(1,1000);
RAAC_new = zeros(1,1000);
RPAC_new = zeros(1,1000);
RCFC = zeros(1,1000);
RAAC = zeros(1,1000);
p_PAC = zeros(1,1000); p_AAC = zeros(1,1000); p_CFC = zeros(1,1000);
p_PAC_new = zeros(1,1000); p_AAC_new = zeros(1,1000);


for j = 1:1000
    [XX] = simfun(mval,mval,'pink','empirical','none',.05);
    RPAC(j) = XX.rpac;
    RPAC_new(j) = XX.rpac_new;
    RCFC(j) = XX.rcfc;
    RAAC(j) = XX.raac;
    RAAC_new(j) = XX.raac_new;
    p_PAC(j) = P.rpac; p_AAC(j) = P.raac; p_CFC(j) = P.rcfc;
    p_PAC_new(j) = P.rpac_new; p_AAC_new(j) = P.raac_new;
end

strname = ['CFC_Simulations_',id];
save(strname)