function [val]=sk_gaussian(N,G,LAM,FA,K)

GK11=real(gkgauss(N,G,LAM,FA,1,1,K));
GK12=real(gkgauss(N,G,LAM,FA,1,2,K));
GK22=real(gkgauss(N,G,LAM,FA,2,2,K));

DET=GK11.*GK22-power(GK12,2);

GK11INV=GK22./DET;
GK22INV=GK11./DET;
GK12INV=-GK12./DET;

val=GK11INV+GK22INV-2*GK12INV;