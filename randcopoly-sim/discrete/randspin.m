clear
close all

K0=1e-2;
KF=1e2;
NK=100;
K=transpose(logspace(log10(K0),log10(KF),NK));


N=20;
G=10;
LAM=-0.1;
FA=0.5;
EPS=10;

CHI=0./(4*FA*(1-FA)*G);

ORD=20;         % (Number of poles-1) in Laplace inversion
ResLayer=100;   % numeber of layers used in the continued fraction form
ORDmax=20;

GK11=real(gk(N,G,LAM,FA,EPS,1,1,K,ORDmax,ORD,ResLayer));
GK12=real(gk(N,G,LAM,FA,EPS,1,2,K,ORDmax,ORD,ResLayer));
GK22=real(gk(N,G,LAM,FA,EPS,2,2,K,ORDmax,ORD,ResLayer));

DET=GK11.*GK22-power(GK12,2);

GK11INV=GK22./DET;
GK22INV=GK11./DET;
GK12INV=-GK12./DET;

GAM=GK11INV+GK22INV-2*GK12INV-2*CHI;

GAM=power(GAM,-1);
