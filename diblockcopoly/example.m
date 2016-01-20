clear;close all;

% Figure 1: make a mean-field phase diagram
N=1e6;
FAV=linspace(0.1,0.499,41);
plotphase(N,FAV)

% Figure 2: make a phase diagram with density fluctuations
N=1e9;
Nbar=1e6;
FAV=linspace(0.3,0.499,41);
plotphaseRG(N,Nbar,FAV);

% % Figure 3: density-density correlations
% N=1e5;FA=0.25;k=logspace(-2,2,50)/sqrt(r2(N));
% G=gamma2(N,FA,k,0);
% figure;plot(k*sqrt(R2(N)),1./G);
% xlabel('(kR)^2');ylabel('1/\Gamma_2')

% Figure 4: vertex functions
N=1e5;
FAV = linspace(0.2,0.5,20);
GAM3 = zeros(length(FAV),1);
GAM4 = zeros(length(FAV),1);
ifa = 1;
for FA = FAV
    [~,ks]=spinodal(N,FA);
    GAM3(ifa)=gamma3(N,FA,ks);
    
    Q1(1:3)=[1,0,0];
    Q2(1:3)=transpose(rotz(pi)*Q1(1:3)');
    Q3(1:3)=-Q2;
    Q4(1:3)=-Q1;
    GAM4(ifa)=gamma4(N,FA,ks,Q1,Q2,Q3,Q4);
    ifa = ifa+1;
end
figure;plot(FAV,-GAM3*N);xlim([0.2,0.5])
figure;plot(FAV,GAM4*N);xlim([0.3,0.5])
filename = sprintf('data/N1e%dgamma.mat',log10(N));
save(filename,'N','GAM3','GAM4','FAV')
