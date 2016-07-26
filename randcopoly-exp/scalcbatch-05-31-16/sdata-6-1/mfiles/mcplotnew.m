clear;

figure;hold
set(gca,'fontsize',18)

%%%%%%%%% INPUT %%%%%%%%%%
testnum=1;chemcp=2;ploton=0;
%chilistname=sprintf('rand-pt-01-08-15-%d/',testnum);
chilistname=sprintf('rand-pt-01-22-15-%d-%d/',testnum,chemcp);
SPINODAL=load('/tower3/home/shifan/wlc-phase-quad-11-27-14/chivals');
EPS=0.01;LAM=-0.75;FA=0.50;

CHIV=load([chilistname,'chilist']);
CHIV=CHIV(1:5:16);
CHIV=CHIV(1:2);

colcnt=1;
NCHI=length(CHIV);
%%%%%%%%% INPUT END %%%%%%%%%%

%%%% start plotting %%%%
ind=find(SPINODAL(:,1)==EPS & SPINODAL(:,2)==LAM); 
CHIS = SPINODAL(ind,3);

for cnum=1:NCHI
CHI=CHIV(cnum);
col=(cnum-1)/(NCHI-1);filename=sprintf('Sdata/SMC_SIM%dCHEM%dCHI%.8f',testnum,chemcp,CHI);

S=0;
for chemcp=1:4
  filename=sprintf('Sdata/SMC_SIM%dCHEM%dCHI%.8f',testnum,chemcp,CHI);
  S=S+load(filename);
end

S=S./4;

%errorbar(S(:,1),S(:,2),S(:,3),'o','linewidth',1,'color',[col 0 1-col])
plot(S(:,1),S(:,2),'s','linewidth',1,'color',[col 0 1-col])
colcnt=colcnt+1;
end

colcnt=1;
for cnum=1:NCHI
CHI=CHIV(cnum);
col=(cnum-1)/(NCHI-1);
  if CHI<CHIS
    filename=sprintf('data/S_FA%2d_EPS%d_LAM%2d_CHI%.8f',FA*100,EPS*100,LAM*100,CHI)
    Sa=load(['/tower3/home/shifan/wlc-phase-quad-11-27-14/',filename]);
plot(Sa(:,1),Sa(:,2),'--','color',[col 0 1-col],'linewidth',1)
 else
filename=sprintf('data/S_FA%2d_EPS%d_LAM%2d_CHI%.8f',FA*100,EPS*100,LAM*100,CHIS)
end
colcnt=colcnt+1;
end


%%%% end plotting %%%%%%%

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('k R_b')
ylabel('S(k)')
%ylim([1e-2,1e4])
ylim([1e-1,1e0])
title(['\chi=',num2str(CHI)])

Ree=2;LBOX=20;DEL=1;
plot([1,1]*2*pi/LBOX*Ree,[1e-1,1e2],'k--','linewidth',2)
plot([1,1]*2*pi/(LBOX/2)*Ree,[1e-1,1e2],'b--','linewidth',2)
plot([1,1]*2*pi/DEL*Ree,[1e-1,1e2],'k--','linewidth',2)
xlim([2*pi/LBOX*Ree/2,2*pi/DEL*Ree])
