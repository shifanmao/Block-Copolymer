%  calculate one-loop density-density correlation
clear;
%close all

u = 1;  % interaction parameter
c = 1;  % chain concentration

N=1e4;
k=logspace(0,5,100)';

gam2 = gamma2(k,N);
loop0 = power(1/u+c*N^2*gam2,-1);
% s = power(1./gam2+c*u*N^2,-1);


figure;
loglog(k,s)
xlabel('k2l_p');ylabel('S(k)')