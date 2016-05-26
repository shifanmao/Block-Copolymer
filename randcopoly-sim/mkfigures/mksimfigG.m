clear;close all
NM=0.05;
R2=-0.5+0.5*exp(-2*NM)+NM;

f1=figure;set(gca,'fontsize',30)
set(f1,'position',[0,0,500,2*600])

f2=figure;hold;set(gca,'fontsize',30)
set(f2,'position',[0,0,800,600])

f3=figure;hold;set(gca,'fontsize',30)
set(f3,'position',[0,0,800,600])

figure(f1);subplot(2,1,1);hold
% plotsim(0.01,-0.75,1,1,1);
[chis,chiv,KS_MF,KS_SIM,SINV_MF,SINV_SIM,D2S_MF,D2S_SIM,ERR]=plotsim(0.01,-0.75,1,1,1);
DEL=1.;plot([4*pi/DEL,4*pi/DEL],[0.05,600],'k--','linewidth',2)
xlim([.6,4*pi/0.5]);ylim([5e-2,200])
set(gca,'fontsize',18)
xlabel('R_Mq');ylabel('<\psi(q)\psi(-q)>')

G=5;
figure(f2);plot(chiv*G,SINV_MF*G,'k-','linewidth',2);
errorbar(chiv*G,SINV_SIM*G,ERR(:,2)*G,'r.-','linewidth',2)

figure(f3);
errorbar(chiv*G,KS_SIM,ERR(:,1),'r.-','linewidth',2)

figure(f1);subplot(2,1,2);hold
% plotsim(0.01,-0.75,1,1,111);
[chis,chiv,KS_MF,KS_SIM,SINV_MF,SINV_SIM,D2S_MF,D2S_SIM,ERR]=plotsim(0.01,-0.75,1,1,111);
DEL=0.8;plot([4*pi/DEL,4*pi/DEL],[0.05,600],'k--','linewidth',2)
xlim([.6,4*pi/0.5]);ylim([5e-2,200])
set(gca,'fontsize',18)
xlabel('R_Mq');ylabel('<\psi(q)\psi(-q)>')

G=10;
figure(f2);
plot(chiv*G,SINV_MF*G,'k-','linewidth',2);
errorbar(chiv*G,SINV_SIM*G,ERR(:,2)*G,'b.-','linewidth',2)

figure(f3);
errorbar(chiv*G,KS_SIM,ERR(:,1),'b.-','linewidth',2)



% clear;close all
% f1=figure;set(gca,'fontsize',30)
% set(f1,'position',[0,0,500,2*600])
% 
% subplot(2,1,1);hold
% plotsim(1,-0.75,1,15,2);DEL=1.;
% plot([4*pi/DEL,4*pi/DEL],[0.05,600],'k--','linewidth',2)
% xlim([.6,4*pi/0.5]);ylim([5e-2,200])
% set(gca,'fontsize',18)
% xlabel('R_Mq');ylabel('<\psi(q)\psi(-q)>/G')
% 
% subplot(2,1,2);hold
% plotsim(1,-0.75,1,15,111);DEL=0.8;
% plot([4*pi/DEL,4*pi/DEL],[0.05,600],'k--','linewidth',2)
% xlim([.6,4*pi/0.5]);ylim([5e-2,200])
% set(gca,'fontsize',18)
% xlabel('R_Mq');ylabel('<\psi(q)\psi(-q)>/G')
% 

% 
% 
% 
% 
% % PLOT ONE DATASET OF STRUCTURE FACTOR
% clear;
% f1=figure;set(gca,'fontsize',30)
% set(f1,'position',[0,0,500,2*600])
% 
% 
% 
% subplot(2,1,1);hold
% set(gca,'fontsize',20)
% vfac = 1;
% lfac = vfac^(1/3);
% 
% % simulation constants
% FA=0.5;  % fraction of A blocks
% M=8;        % number of blocks
% G=5*vfac;  % number of discrete monomers
% EPS=1/vfac;
% LAM=-0.75;
% NM=G*EPS;  % number of Kuhn steps per monomer
% 
% PLOTON=0;
% chiind=1;
% 
% [KS_SIM,SINV_SIM,D2S_SIM,S,CHI]=plotsim5(EPS,LAM,vfac,PLOTON,15,2);
% 
% %%%% limit simulation range of structure factor %%%%
% BOXL=20/lfac;DEL=1/lfac;RM=2;
% SIND1 = find(S(:,1) >= RM*2*pi/BOXL);SIND1 = SIND1(1);
% SIND2 = find(S(:,1) >= RM*2*pi/DEL);SIND2 = SIND2(1);
% SRANGE = SIND1:SIND2;
% S = S(SRANGE,:);
% % S(1,1) = S(2,1)-1e-2;
% % S(end,1) = S(end-1,1)+1e-2;
% % S(1,end) = 0; S(end,end) = 0;
% S = [S(1,1)-1e-2,0 ; S; S(end,1)+1e-2,0];
% 
% % numerical integration of simulation structure factor
% % p = patch(S(:,1),S(:,2),'b','LineWidth',1.5);
% 
% % PLOT MEAN-FIELD STRUCTURE FACTOR
% R2=-0.5+0.5*exp(-2*NM)+NM;
% k=logspace(-1,3,100)./sqrt(R2); % wavevectors
% val = s2invwlc(M,NM,FA,LAM,k);
% plot(k.*sqrt(R2),1./(-2*CHI+EPS*val)/G,'k--','linewidth',3)            
% 
% % PLOT SIMULATION STRUCTURE FACTOR
% plot(S(:,1),S(:,2)/G,'k.-','linewidth',3,'markersize',20)
% 
% % PLOT APPROXIMATED STRUCTURE FACTOR
% x = logspace(-1,3,100);
% y = 1./(SINV_SIM+D2S_SIM/G*(x-KS_SIM).^2);
% plot(x,y/G,'k-','linewidth',2)
% set(gca,'yscale','linear');
% set(gca,'xscale','linear');
% xlim([0,20]);ylim([0,.5]);
% % xlim([2e-1,5e1]);ylim([0,3]);
% 
% xlabel('R_Mq');
% ylabel('<\psi(q)\psi(-q)>/G')
% % xlim([.6,12]);ylim([1e-1,2e2]);
% % set(gca,'Ytick',[1e-1 1e0 1e1 1e2])
% % set(gca,'Xtick',[1,10])
% % set(gca,'xscale','log');
% % legend('Mean-field','Simulation','Fredrickson-Helfand');
% legend('Mean-field','Simulation','Lorentzian Fit');
% box on
% 
% % EVALUATE INTEGRALS AND SUMMATIONS
% xsum = S(:,1)/RM;
% ysum = xsum.^2.*S(:,2);
% SUM_SIM = 4*pi*trapz(xsum,ysum)/G
% INT_SIM = 4*pi^2*(KS_SIM/RM)^2/RM * sqrt(G/SINV_SIM/D2S_SIM)/G
% KS_SIM
% 1/(G*SINV_SIM)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% subplot(2,1,2);hold
% set(gca,'fontsize',20)
% vfac = 2;
% lfac = vfac^(1/3);
% 
% % simulation constants
% FA=0.5;  % fraction of A blocks
% M=8;        % number of blocks
% G=5*vfac;  % number of discrete monomers
% EPS=1/vfac;
% LAM=-0.75;
% NM=G*EPS;  % number of Kuhn steps per monomer
% 
% PLOTON=0;
% chiind=1;
% 
% [KS_SIM,SINV_SIM,D2S_SIM,S,CHI]=plotsim5(EPS,LAM,vfac,PLOTON,15,111);
% 
% %%%% limit simulation range of structure factor %%%%
% BOXL=20/lfac;DEL=1/lfac;RM=2;
% SIND1 = find(S(:,1) >= RM*2*pi/BOXL);SIND1 = SIND1(1);
% SIND2 = find(S(:,1) >= RM*2*pi/DEL);SIND2 = SIND2(1);
% SRANGE = SIND1:SIND2;
% S = S(SRANGE,:);
% % S(1,1) = S(2,1)-1e-2;
% % S(end,1) = S(end-1,1)+1e-2;
% % S(1,end) = 0; S(end,end) = 0;
% S = [S(1,1)-1e-2,0 ; S; S(end,1)+1e-2,0];
% 
% % numerical integration of simulation structure factor
% % p = patch(S(:,1),S(:,2),'b','LineWidth',1.5);
% 
% % PLOT MEAN-FIELD STRUCTURE FACTOR
% R2=-0.5+0.5*exp(-2*NM)+NM;
% k=logspace(-1,3,100)./sqrt(R2); % wavevectors
% val = s2invwlc(M,NM,FA,LAM,k);
% plot(k.*sqrt(R2),1./(-2*CHI+EPS*val)/G,'k--','linewidth',3)            
% 
% % PLOT SIMULATION STRUCTURE FACTOR
% plot(S(:,1),S(:,2)/G,'k.-','linewidth',3,'markersize',20)
% 
% % PLOT APPROXIMATED STRUCTURE FACTOR
% x = logspace(-1,3,100);
% y = 1./(SINV_SIM+D2S_SIM/G*(x-KS_SIM).^2);
% plot(x,y/G,'k-','linewidth',2)
% set(gca,'yscale','linear');
% set(gca,'xscale','linear');
% xlim([0,20]);ylim([0,.5]);
% % xlim([2e-1,5e1]);ylim([0,3]);
% 
% xlabel('R_Mq');
% ylabel('<\psi(q)\psi(-q)>/G')
% % xlim([.6,12]);ylim([1e-1,2e2]);
% % set(gca,'Ytick',[1e-1 1e0 1e1 1e2])
% % set(gca,'Xtick',[1,10])
% % set(gca,'xscale','log');
% % legend('Mean-field','Simulation','Fredrickson-Helfand');
% legend('Mean-field','Simulation','Lorentzian Fit');
% box on
% 
% % EVALUATE INTEGRALS AND SUMMATIONS
% xsum = S(:,1)/RM;
% ysum = xsum.^2.*S(:,2);
% SUM_SIM = 4*pi*trapz(xsum,ysum)/G
% INT_SIM = 4*pi^2*(KS_SIM/RM)^2/RM * sqrt(G/SINV_SIM/D2S_SIM)/G
% KS_SIM
% 1/(G*SINV_SIM)