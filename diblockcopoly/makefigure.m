%% Makes the plots of 'Semiflexible diblock' paper
clear;close all;

% Figure 1: density-density correlations
N=1e5;FA=0.25;k=logspace(-2,2,50)/sqrt(r2(N));
G=gamma2(N,FA,k,0);
figure;plot(k*sqrt(R2(N)),1./G);
xlabel('(kR)^2');ylabel('1/\Gamma_2')

% Figure 2: vertex functions
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

% % Figure 3: spinodal
% FAV = 0.5;
% NV = logspace(-1,3,21);
% ks = zeros(length(NV),1);
% chis = zeros(length(NV),1);
% inm = 1;
% for N = NV
%     [chis(inm),ks(inm)]=spinodal(N,FAV);
%     inm = inm+1;
% end
% figure;semilogx(NV,chis.*NV')
% figure;loglog(NV,1./ks)
% 
% Figure 4: phase diagram
[chis,chi13,chi36,chi12,chi23,chi26]=plotphase(N,linspace(0.15,0.5,20));