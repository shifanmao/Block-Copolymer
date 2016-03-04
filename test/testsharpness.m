close all;clear

N=100;
NM=100;
FA=0.5;
LAM=-0.99;
CHI=0;

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-2;  % minimum wavevector
KF=1e1;   % maximum wavevector
NK=2001;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

S2INV=s2invgc(N,NM,FA,LAM,K);
figure;hold;set(gca,'fontsize',20)
loglog(RM*K,1./(-2*CHI+S2INV),'k-','linewidth',2);
xlabel('R_Mq');ylabel('S(q)')

SMAX=max(1./S2INV);
INDMAX=find((SMAX==1./S2INV));
KMAX=RM*K(INDMAX);
set(gca,'yscale','linear');set(gca,'xscale','linear')

D2S = 0.013409;
x=RM*K(INDMAX-100:INDMAX+100);
y=SMAX-0.5*D2S*NM*(x-KMAX).^2;
plot(x,y,'k--','linewidth',2)
plot(KMAX,SMAX,'k.','markersize',20)
legend('<\psi^2>','local fit')