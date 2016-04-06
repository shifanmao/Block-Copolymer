clear;
addpath('functions/')
addpath('chainstats/')
addpath('misc/')

N=100;
NM=1e1;
LAM=-0.75;
FA= 0.5;
CHI=0;
k=logspace(-1,3,100)/sqrt(r2(NM));

% Mean-field density-density correlation
Gmf=gamma2(N,NM,LAM,FA,k,CHI);
Smf=1./(Gmf);  % Density-density correlation
figure;loglog(k*sqrt(r2(NM)),Smf)

% Spinodal
N=100;
FAV=linspace(0.1,0.9,81);  % range of A monomer chemical composition
NMV=1e4;
LAM=0.;
[chis,ks,d2gam2]=spinodal(N,NMV,LAM,FAV);
figure;plot(FAV,chis.*NMV);

% Vertex functions
N=100;
NM=1e5;
NQ=1;  % number of wavevector sets in calculating GAM4
LAM=0.;
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
[gam3,gam4]=calcgamma(N,NM,LAM,FAV,NQ);

figure;plot(FAV,-gam3*NM,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)')
figure;plot(FAV,gam4*NM,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)')