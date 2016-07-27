% calculates zero-loop density-density correlation
clear;
%close all

u = 1;  % interaction parameter
c = 1;  % chain concentration

N=1e0;
k=logspace(-2,5,100)';

gam2 = gamma2(k,N);
loop0 = power(1/u+c*N^2*gam2,-1);
s = u^(-1)-u^(-2)*loop0;

figure;
loglog(k,s)
xlabel('k2l_p');ylabel('S(k)')