%MAKE VLO ARRAY time by numtrials
%MAKE VHI ARRAY time by numtrials
%compute MI between random column of VLO, random column of VHI
%sample with replacement
%200 surrogate MI values

L = 100;
Vlo_array = []; Vhi_array = [];

for i = 1:L
    [XX,P,Vlo,Vhi,t] = simfun(1,0,'GLM','theoretical',.05);
    Vlo_array = [Vlo_array,Vlo'];
    Vhi_array = [Vhi_array,Vhi'];
end

num_samples = 200;

k1 = ceil(rand(1,num_samples)*L);
k2 = ceil(rand(1,num_samples)*L);

MI = zeros(1,num_samples);
for i = 1:num_samples
    MI(i) = modulation_index(Vlo_array(:,k1(i)),Vhi_array(:,k2(i)));
end

% 4.6809e-04
