clear;

NM_WLC=1e1;
QM=logspace(-1.5,3.5,100)';
ang = pi/3;
gam3 = gamma3(QM,ang,NM_WLC);

% % load option
% NM_WLC=1e1;
% filename = sprintf('gam3N1e%d',log10(NM_WLC));
% x = load(filename);
% QM = 10.^(x(:,1));
% gam3 = 10.^(x(:,2));

figure('Position', [100, 100, 1200, 900])
hold;set(gca,'fontsize',50)
plot(QM,gam3,'b-','linewidth',6);
plot(QM,4./power(QM,4),'k--','linewidth',4);
xlim([10^-1.5,10^2]);ylim([1e-6,1]);box on 
set(gca,'xscale','log');set(gca,'yscale','log');
xlabel('Wavevector kR_g');ylabel('Cubic-order Vertex \Gamma^{(3)}')

% save to figure
filename = sprintf('gam3N1e%d',log10(NM_WLC));
saveas(gcf,strcat(filename,'.eps'),'epsc')

% save to data
dlmwrite(filename,[log10(QM*Rg),log10(gam3)],'precision','%.2f')