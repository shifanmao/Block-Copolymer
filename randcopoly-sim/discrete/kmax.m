function [kval,sval]=kmax(N,G,LAM,FA,EPS)
ORDmax=30;
ORD=30;
ResLayer=100;

DK=1e-1;

NG=EPS*G;
R2=r2wlc(NG);

[kval,sval] = fminbnd(@(K) sk(N,G,LAM,FA,EPS,K,ORDmax,ORD,ResLayer), -0.5/sqrt(R2), 10/sqrt(R2), optimset('TolX',1e-8,'Display','off'));

if kval<(1e-3/sqrt(R2))
  kval=0;
  sval= sk(N,G,LAM,FA,EPS,1e-3/sqrt(R2),ORDmax,ORD,ResLayer);
end
