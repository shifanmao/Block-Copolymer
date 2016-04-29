% clear;close all
addpath('../misc/')
addpath('../functions/')
addpath('../chainstats/')
addpath('../chainstats/integrals/')

N=1e3;
NM=100;
C=1e2;
LAM=-0.75;
FA=0.5;
densityRG(N,C,NM,LAM,FA)
% 
% LAMV=-0.9:0.10:0.9;
% gam4=zeros(length(LAMV),1);
% gam4rep=zeros(length(LAMV),1);
% NQ=1;
% cnt=1;
% for LAM=LAMV
%      [~,gam4(cnt),gam4rep(cnt)]=calcgamma(N,NM,LAM,FA,NQ);
%     cnt=cnt+1;
% end

LAMV=-0.9:0.01:0.9;
chit=zeros(length(LAMV),1);
chis=zeros(length(LAMV),1);
phase=zeros(length(LAMV),1);
cnt=1;
for LAM=LAMV
    [chit(cnt),phase(cnt)]=spinodalRG(N,C,NM,LAM,FA);
    [chis(cnt),ks,d2gam2]=spinodal(N,NM,LAM,FA);
    cnt=cnt+1;
end