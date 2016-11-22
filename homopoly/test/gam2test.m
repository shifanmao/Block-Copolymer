clear;

NM_WLC=40;
Rg = sqrt(NM_WLC/6);
k=logspace(-1,2,100)';
gam2 = gamma2(k,NM_WLC);

s2g = 2./(k*Rg).^4.*(exp(-(k*Rg).^2)+(k*Rg).^2-1);
s2gasym = 2./(k*Rg).^2;
s2wlcasym = .54./(k*Rg^2).^1;
s2expand = 1-(k*Rg).^2/3;

figure('Position', [100, 100, 1200, 900])
hold;set(gca,'fontsize',50)

% plot simulation results
N = 40;CHI = 0.0;
filename = sprintf('N%.2fCHI%.2f',N,CHI);
t = load(strcat('../data/MCsim/',filename));
plot(t(:,1),t(:,2),'^','MarkerSize',20,'linewidth',2)

% plot analytical theory
plot(k,gam2,'b-','linewidth',6);
plot(k,s2g,'k--','linewidth',4)
plot(k,s2wlcasym,'k--','linewidth',4)
% plot(k,s2expand,'r--','linewidth',4)
xlim([1e-1,10^1.5]);ylim([10^-2.5,1e0]);box on
set(gca,'xscale','log');set(gca,'yscale','log')
title(strcat('u_0=',sprintf('%.1f',CHI)))
xlabel('Wavevector K');ylabel('Quadratic-order Vertex \Gamma^{(2)}')

% save to figure
%filename = sprintf('gam2N%d',NM_WLC);
saveas(gcf,strcat(filename,'.eps'),'epsc')