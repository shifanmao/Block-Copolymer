% clear;
% 
% N=1e4;NM=1e2;FA=0.5;
% % kv=logspace(-1,5,20)/N/NM;
% k=1e-2/N/NM;
% LAMV=0;
% gam4rep=zeros(length(LAMV),1);
% gam4   =zeros(length(LAMV),1);
% Q1=[1,0,0];Q2=-Q1;Q3=Q1;Q4=-Q1;
% 
% cnt=1;
% for LAM=LAMV
%     LAM
%     gam4(cnt)=gamma4(N,NM,LAM,FA,k,Q1,Q2,Q3,Q4);
%     gam4rep(cnt)=gamma4rep(N,NM,LAM,FA,k,Q1,Q2,Q3,Q4);
%     cnt=cnt+1;
% end
% 
% figure;hold
% semilogx(k*N*NM,gam4,'k-')
% semilogx(k*N*NM,gam4rep,'k--')


% clear;
% LAM=1;
% NQ=1;
% N=100;NM=1e4;
% FAV=linspace(0.1,0.5,20);
% [gam3,gam4,gam4rep]=calcgamma(N,NM,LAM,FAV,NQ);
% figure;plot(FAV,gam4*NM,FAV,gam4rep*NM)

% 
% clear;
% NQ=1;
% N=1e2;NM=1e2;
% FA=0.5;
% 
% LAMV=-0.99:0.01:0.99;
% gam3   =zeros(length(LAMV),1);
% gam4   =zeros(length(LAMV),1);
% gam4rep=zeros(length(LAMV),1);
% cnt=1;
% for LAM=LAMV
%     LAM
%     [gam3,gam4(cnt),gam4rep(cnt)]=calcgamma(N,NM,LAM,FA,NQ);
%     cnt=cnt+1;
% end
% 
% figure;plot(LAMV,-gam3*NM)
% figure;plot(LAMV,gam4*NM,LAMV,gam4rep*NM)
% legend('gam4','gam4rep')

% clear;
% LAM=.99;FA=0.5;
% N=1e4;NM=1e2;
% k=1e-4;Q1=[1,0,0]*1e-6;Q2=-Q1;Q3=Q1;Q4=-Q1;
% gamma4rep(N,NM,LAM,FA,k,Q1,Q2,Q3,Q4)*N*NM

clear;
N=1e3;NM=1e2;
FA=0.5;
LAM=0;

cnt=1;
kv=logspace(0,5,20);
for k=kv/N/NM
    k
    Q1=[1,0,0]*k;Q2=-Q1;Q3=Q1;Q4=-Q1;
    S4=s4gcrep(N,NM,LAM,FA,Q1,Q2,Q3,Q4);
    S2=s2gc(N,NM,LAM,FA,norm(Q1));
    
    s4aaaa(cnt)=S4(1,1,1,1)/power(N*NM,4);
    s22aa(cnt)=S2(1,1)/power(N*NM,2);
    cnt=cnt+1;
end

figure;plot(kv,s22aa,kv,s4aaaa)
set(gca,'xscale','log');set(gca,'yscale','linear');