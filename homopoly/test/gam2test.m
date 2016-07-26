clear;

NM_WLC=1e2;Rg = sqrt(NM_WLC/6);
k=logspace(0,2,100)'/Rg;
gam2 = gamma2(k,NM_WLC);

s2g = 2./(k*Rg).^4.*(exp(-(k*Rg).^2)+(k*Rg).^2-1);
s2gasym = 2./(k*Rg).^2;
s2wlcasym = .54./(k*Rg^2).^1;

figure('Position', [100, 100, 1200, 900])
hold;set(gca,'fontsize',50)
plot(k*Rg,gam2,'b-','linewidth',6);
plot(k*Rg,s2g,'k--','linewidth',4)
% plot(k*Rg,s2gasym,'r--')
plot(k*Rg,s2wlcasym,'k--','linewidth',4)
xlim([1e0,1e2]);ylim([1e-3,1e0]);box on
set(gca,'xscale','log');set(gca,'yscale','log')
xlabel('Wavevector kR_g');ylabel('Quadratic-order Vertex \Gamma^{(2)}')

filename = sprintf('gam2N1e%d',log10(NM_WLC));

% save to figure
saveas(gcf,strcat(filename,'.eps'),'epsc')