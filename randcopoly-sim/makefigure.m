%% Makes the plots of 'Impact of' paper
clear;close all;

%% Pre-calculation
% find lifshitz point at different NM
load('data/spinwlc.mat');
LAML = zeros(length(NMV),1);
for ii = 1:length(NMV)
    ind = find(ks(:,ii)<=0);
    LAML(ii)=LAMV(ind(1));
end

%% Figure 1: density-density correlations
NM=10;k=logspace(-2,2,50)/sqrt(r2(NM));
G=gamma2(100,NM,-.75,0.5,k,0.);
figure;loglog(k*sqrt(r2(NM)),1./G);

%% Figure 2:  spinodal
figure;hold

load('data/spingc.mat');
plot(LAMV,chis*NMV)

load('data/spinwlc.mat');
for ii = 1:length(NMV)
    plot(LAMV,chis(:,ii)*NMV(ii))
end

%% Figure 3: critical wavemode
figure;hold

load('data/spingc.mat');
plot(LAMV,ks*sqrt(r2(NMV)))

load('data/spinwlc.mat');
for ii = 1:length(NMV)
    plot(LAMV,ks(:,ii).*sqrt(r2(NMV(ii))))
end

%% Figure 4: peak sharpness
figure;hold
load('data/spingc.mat');
plot(LAMV-LAML(end),log(abs(d2gam2)))

load('data/spinwlc.mat');
for ii = 1:length(NMV)
    plot(LAMV-LAML(ii),log(abs(d2gam2(:,ii))))
end

%% Figure 5: Sinv(q*)
LAM=-0.75;
EPSV=[0.01,0.10,1.00];
color = ['r','b','k'];

f2=figure;hold;set(gca,'fontsize',30)
set(f2,'position',[10,0,800,600])
cnt = 1;
for EPS=EPSV
    [chis,chiplot,ks,ksim,sinv_theory,sinv_sim]=plotmcsim(EPS,LAM,0);
    plot(chiplot,sinv_sim,'o','color',color(cnt),'linewidth',3)
    cnt = cnt+1;
end

cnt = 1;
for EPS=EPSV
    [chis,chiplot,ks,ksim,sinv_theory,sinv_sim]=plotmcsim(EPS,LAM,0);
    plot(chiplot,sinv_theory,'-','color',color(cnt),'linewidth',3)
    cnt = cnt+1;
end

xlabel('\chivG');ylabel('S^{-1}(q^*)')
ylim([0,2]);xlim([0,20]);box on
legend('N_M=0.05','N_M=0.50','N_M=5.00')

%% Figure 6: q*
f2=figure;hold;set(gca,'fontsize',30)
set(f2,'position',[10,0,800,600])
cnt = 1;
for EPS=EPSV
    [chis,chiplot,ks,ksim,sinv_theory,sinv_sim]=plotmcsim(EPS,LAM,0);
    plot(chiplot,ksim,'o-','color',color(cnt),'linewidth',3)
    %chis*5*EPS
    cnt = cnt+1;
end

cnt = 1;
for EPS=EPSV
    [chis,chiplot,ks,ksim,sinv_theory,sinv_sim]=plotmcsim(EPS,LAM,0);
    plot(chiplot,repmat(ks*sqrt(r2(EPS*5)),length(chiplot),1),'-','color',color(cnt),'linewidth',3)
    %chis*5*EPS
    cnt = cnt+1;
end

xlabel('\chivG');ylabel('q^* R_M')
%ylim([0,2]);
xlim([0,20]);
legend('\epsilon=0.01','\epsilon=0.10','\epsilon=1.00')