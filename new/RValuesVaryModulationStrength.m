% AAC
MOD = [0:.05:.95];
RPAC = zeros(20,1000); PPAC = zeros(20,1000);
RCFC = zeros(20,1000); PCFC = zeros(20,1000);
RAAC = zeros(20,1000); PAAC = zeros(20,1000);

for i = 1:20
    mval = MOD(i);
    for j = 1:1000
        [XX,P] = simfun_rsquared(0,mval,'pink','empirical','none',.05);
        RPAC(i,j) = XX.rpac; PPAC(i,j) = P.rpac;
        RCFC(i,j) = XX.rcfc; PCFC(i,j) = P.rcfc;
        RAAC(i,j) = XX.raac; PAAC(i,j) = P.raac;
    end
end
strname = ['Emp_AAC_Vary_Mod'];
save(strname)

% PAC
MOD = [0:.05:.95];
RPAC = zeros(20,1000); PPAC = zeros(20,1000);
RCFC = zeros(20,1000); PCFC = zeros(20,1000);
RAAC = zeros(20,1000); PAAC = zeros(20,1000);

for i = 1:20
    mval = MOD(i);
    for j = 1:1000
        [XX,P] = simfun_rsquared(mval,0,'pink','empirical','none',.05);
        RPAC(i,j) = XX.rpac; PPAC(i,j) = P.rpac;
        RCFC(i,j) = XX.rcfc; PCFC(i,j) = P.rcfc;
        RAAC(i,j) = XX.raac; PAAC(i,j) = P.raac;
    end
end
strname = ['Emp_PAC_Vary_Mod'];
save(strname)

% PAC and AAC
MOD = [0:.05:.95];
RPAC = zeros(20,1000); PPAC = zeros(20,1000);
RCFC = zeros(20,1000); PCFC = zeros(20,1000);
RAAC = zeros(20,1000); PAAC = zeros(20,1000);

for i = 1:20
    mval = MOD(i);
    for j = 1:1000
        [XX,P] = simfun_rsquared(mval,mval,'pink','empirical','none',.05);
        RPAC(i,j) = XX.rpac; PPAC(i,j) = P.rpac;
        RCFC(i,j) = XX.rcfc; PCFC(i,j) = P.rcfc;
        RAAC(i,j) = XX.raac; PAAC(i,j) = P.raac;
    end
end
strname = ['Emp_PAC_AAC_Vary_Mod'];
save(strname)