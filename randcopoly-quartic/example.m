clear;
addpath('functions/')
addpath(genpath('chainstats/'))
addpath('misc/')

N=100;
NM=1e3;
LAMV=linspace(-1,1,41);
FA=0.5;
CHI=0;
k=logspace(-1,3,100)/sqrt(r2(NM));

chis = zeros(length(LAMV),1);
ks = zeros(length(LAMV),1);
gam3 = zeros(length(LAMV),1);
gam4 = zeros(length(LAMV),1);
gam4rep = zeros(length(LAMV),1);
cnt = 1;
for LAM = LAMV
    % Spinodal
    [chis(cnt),ks(cnt),d2gam2]=spinodal(N,NM,LAM,FA);

    % Vertex functions
    NQ=1;  % number of wavevector sets in calculating GAM4
    [gam3(cnt),gam4(cnt),gam4rep(cnt)]=calcgamma(N,NM,LAM,FA,NQ);

    cnt = cnt+1;
end
figure;plot(LAMV,chis*NM);

figure;plot(LAMV,-gam3*NM,'k-','linewidth',2);
xlabel('\lambda');ylabel('-N\Gamma_3(q^*)')
figure;plot(LAMV,gam4*NM,'k-','linewidth',2);
xlabel('\lambda');ylabel('N\Gamma_4(q^*)')
figure;plot(LAMV,-gam4rep*NM,'kx-','linewidth',2);
xlabel('\lambda');ylabel('N\Gammarep_4(q^*)')

figure;plot(LAMV,(gam4-gam4rep)*NM,'k-','linewidth',2);
xlabel('f_A');ylabel('N\Gammarep_4(q^*)')




N=2;
NM=1e3;
LAM=-1;
FA= 0.5;
CHI=0;
k=logspace(-1,3,100)/sqrt(r2(NM));

% Mean-field density-density correlation
Gmf=gamma2(N,NM,LAM,FA,k,CHI);
Smf=1./(Gmf);  % Density-density correlation

% Spinodal
FAV=linspace(0.1,0.9,21);  % range of A monomer chemical composition
[chis,ks,d2gam2]=spinodal(N,NM,LAM,FAV);
figure;loglog(k*sqrt(r2(NM)),Smf)
figure;plot(FAV,chis*NM);

% Vertex functions
NQ=1;  % number of wavevector sets in calculating GAM4
[gam3,gam4,gam4rep]=calcgamma(N,NM,LAM,FAV,NQ);

figure;plot(FAV,-gam3*NM,'k-','linewidth',2);xlim([0.2,0.5]);
xlabel('f_A');ylabel('-N\Gamma_3(q^*)')
figure;plot(FAV,gam4*NM,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gamma_4(q^*)')
figure;plot(FAV,gam4rep*NM,'k-','linewidth',2);xlim([0.3,0.5]);
xlabel('f_A');ylabel('N\Gammarep_4(q^*)')