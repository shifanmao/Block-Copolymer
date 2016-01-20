% make a phase diagram with density fluctuations
N=1e9;
Nbar=1e6;
FAV=linspace(0.3,0.499,41);

% calculate renormalized ODT and OOTs
% [chit,chi13,chi36,chi1,chi3,chi6]=plotphaseRG(N,Nbar,FAV);
[chit,phase,chi13,chi36]=plotphaseRG(N,Nbar,FAV);
[chis,~,~]=spinodal(N,FAV);

% plot phase diagram
% ind13=find(chit<chi1);
% ind36=find(chit==chi6);
ind13=find(phase==3 | phase==6);
ind36=find(phase==6);

figure;hold;set(gca,'fontsize',18)
plot(FAV,chit*N,'k','linewidth',1.5)
plot(FAV,chis*N,'k--','linewidth',1.5)
plot(FAV(ind13),chi13(ind13)*N,'r','linewidth',1.2)
plot(FAV(ind36),chi36(ind36)*N,'b','linewidth',1.2)
xlabel('f');ylabel('\chi N')
xlim([FAV(1),FAV(end)]);ylim([10.,15])