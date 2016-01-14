%% Makes the plots of 'Semiflexible diblock' paper
clear;close all;

addpath('chainstats/')
addpath('misc/')

% Figure 1: density-density correlations
N=1e5;
NM=1e5;
FA=0.25;
k=logspace(-2,2,50)/sqrt(r2(N));
G=gamma2(N,NM,FA,k,0);
figure;plot(k*sqrt(r2(N)),1./G);
xlabel('(kR)^2');ylabel('1/\Gamma_2')

% Figure 2: vertex functions
N=1e5;
NM=1e5;
FAV=linspace(0.2,0.5,20);
GAM3=zeros(length(FAV),1);
GAM4=zeros(length(FAV),1);
ifa=1;
for FA=FAV
    [~,ks]=spinodal(N,NM,FA);
    GAM3(ifa)=gamma3(N,NM,FA,ks);

    Q1(1:3)=[1,0,0];
    Q2(1:3)=transpose(rotz(pi)*Q1(1:3)');
    Q3(1:3)=-Q2;
    Q4(1:3)=-Q1;
    GAM4(ifa)=gamma4(N,NM,FA,ks,Q1,Q2,Q3,Q4);
    ifa = ifa+1;
end
figure;plot(FAV,-GAM3*N);xlim([0.2,0.5])
figure;plot(FAV,GAM4*N);xlim([0.3,0.5])

% Figure 3: spinodal
FAV = 0.5;
N = 1e3;
NMV = logspace(-1,3,21);
ks = zeros(length(NMV),1);
chis = zeros(length(NMV),1);
inm = 1;
for NM = NMV
    [chis(inm),ks(inm)]=spinodal(N,NM,FAV);
    inm = inm+1;
end
figure;semilogx(NMV,chis.*N)
figure;loglog(NMV,1./ks)
% 
% % Figure 4: phase diagram
% % [chis,chi13,chi36,chi12,chi23,chi26]=plotphase(N,linspace(0.15,0.5,20));