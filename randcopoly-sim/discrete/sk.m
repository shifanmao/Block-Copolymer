function [val]=sk(N,G,LAM,FA,EPS,K)

ORDmax=30;
ORD=30;
ResLayer=100;

GK11=real(gk(N,G,LAM,FA,EPS,1,1,K,ORDmax,ORD,ResLayer));
GK12=real(gk(N,G,LAM,FA,EPS,1,2,K,ORDmax,ORD,ResLayer));
GK22=real(gk(N,G,LAM,FA,EPS,2,2,K,ORDmax,ORD,ResLayer));

DET=GK11.*GK22-power(GK12,2);

GK11INV=GK22./DET;
GK22INV=GK11./DET;
GK12INV=-GK12./DET;

val=GK11INV+GK22INV-2*GK12INV;