%% calculate structure factors
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')
clear;

FA = 0.5;
k = logspace(-1,2,100);

N = 100;
G=gamma2(N,FA,k,0);

filename = sprintf('SDIB_N%.2f',log10(N));
dlmwrite(filename,[k',1./G])