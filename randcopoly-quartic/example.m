clear;
addpath('functions/')
addpath('chainstats/')
addpath('misc/')

N=200;
NM=1e4;
LAM=0;
FA= 0.5;
CHI=0;
k=logspace(-1,3,100)/sqrt(r2(NM));

% Mean-field density-density correlation
Gmf=gamma2(N,NM,LAM,FA,k,CHI);
Smf=1./(Gmf);  % Density-density correlation

% Spinodal
FAV=linspace(0.1,0.9,81);  % range of A monomer chemical composition
[chis,ks,d2gam2]=spinodal(N,NM,LAM,FAV);
% figure;loglog(k*sqrt(r2(NM)),Smf)
% figure;plot(FAV,chis*NM);

% Vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
FAV=linspace(0.1,0.5,21);  % range of A monomer chemical composition
[gam3,gam4,gam4rep]=calcgamma(N,NM,LAM,FAV,NQ);

figure;plot(FAV,-gam3*NM,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)')
figure;plot(FAV,gam4*NM,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)')
figure;plot(FAV,gam4rep*NM,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gammarep_4(q^*)')