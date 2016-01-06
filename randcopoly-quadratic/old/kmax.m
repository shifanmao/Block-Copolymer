function [kval,sval]=kmax(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer)

DK=0.01;

R2=r2wlc(NM);
%K0=3.3*sqrt(6)/sqrt(R2);
K0=0.5*3.3*sqrt(6)/sqrt(R2);


K=fzero(@(K) (s2invwlc(N,NM,FA,LAM,K+DK,d,ORDmax,ORD,ResLayer)-...
    s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer))/DK,K0);

if sqrt(R2)*K <= 0.01
    kval=0;
else
    A=s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer);
    B=s2invwlc(N,NM,FA,LAM,K+DK/2,d,ORDmax,ORD,ResLayer);
    C=s2invwlc(N,NM,FA,LAM,K+DK,d,ORDmax,ORD,ResLayer);
    K=(A+C-2*B)/(DK/2)^2;
    kval=(A-C)/(K*DK)+K+DK/2;
end

sval=s2invwlc(N,NM,FA,LAM,kval,d,ORDmax,ORD,ResLayer);
