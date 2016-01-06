function [kval,sval]=kmaxwlc(N,NM,FA,LAM,d,ORDmax,ORD,ResLayer)

R2=r2wlc(NM);

NK=50;

KV=transpose(logspace(-2,2,NK))/sqrt(R2);

S=s2invwlc(N,NM,FA,LAM,KV,d,ORDmax,ORD,ResLayer);
[~,IND]=min(S);

if IND ==1
    kval=1e-2/sqrt(R2);
    sval=s2invwlc(N,NM,FA,LAM,kval,d,ORDmax,ORD,ResLayer);
    kval=0;
else
    
    KV2=transpose(linspace(KV(IND-1),KV(IND+1),NK));
    DK=KV2(2)-KV2(1);
    S=s2invwlc(N,NM,FA,LAM,KV2,d,ORDmax,ORD,ResLayer);
    [~,IND]=min(S);
    K=KV2(IND);
    
    A=s2invwlc(N,NM,FA,LAM,K-DK,d,ORDmax,ORD,ResLayer);
    B=s2invwlc(N,NM,FA,LAM,K,d,ORDmax,ORD,ResLayer);
    C=s2invwlc(N,NM,FA,LAM,K+DK,d,ORDmax,ORD,ResLayer);
    KAP=(A+C-2*B)/DK^2;
    kval=(A-C)/(2*KAP*DK)+K;
    
    sval=s2invwlc(N,NM,FA,LAM,kval,d,ORDmax,ORD,ResLayer);
    
end


