%% calculate mean-field spinodal and critical wavelength at FA=0.5
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

clear;
% Figure 3: mean-field spinodal and critical wavelength at FA=0.5
NV=logspace(0,3,16);  % number of statistical steps of total chain
[chis,ks,d2gam2]=spinodal(NV,0.5);

data = [log10(NV)',(chis.*NV)',ks'];
dlmwrite('chivals',data,'precision','%.3f')