clear;

N = 100;
% NM = 1e-2;
NM = 1e2;
PHIP = 0.5;

RM = sqrt(r2(NM));    
KV = logspace(-2, 2, 101)/RM;

FAV = [0.5, 0.8];
LAMV = [-0.75, 0.00];
CHIABV = [0,0.2,0.4,0.6,0.8];

for LAM = LAMV
    for FA = FAV

    figure;hold
    set(gca,'fontsize',18)
    for ii = 1:length(CHIABV)
        COL = (ii-1)/(length(CHIABV)-1)
        CHI = CHIABV(ii);

        [CHIABS, KS, ~, ~] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
        [EIG,EIGV] = gamma2_solvent(N,NM,LAM,FA,KV,[0, CHI; CHI, 0]*CHIABS,PHIP);

        plot(KV*RM, 1./EIG(:,1)/NM, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
        plot(KV*RM, 1./EIG(:,2)/NM, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
    end

    box on
    xlabel('Rq');
    ylabel('$\langle {\psi}^{\bf{\top}}(\vec{q}) {\psi}(-\vec{q}) \rangle/N_{m}v$',...
        'Interpreter', 'Latex')
    xlim([min(KV), max(KV)]*RM)
    set(gca,'fontsize',18)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(strcat('NM=', num2str(NM), ...
                 ' / LAM=', num2str(LAM), ...
                 ' / FA = ', num2str(FA)))
    end
end








% clear;
% 
% N = 100;
% NM = 0.01;
% PHIP = 0.5;
% 
% RM = sqrt(r2(NM));    
% KV = logspace(-2, 2, 51)/RM;
% 
% FA = 0.1;
% LAM = -0.38;
% % CHIABV = [0,0.2,0.4,0.6,0.8];
% CHIABV = linspace(0, 1.01, 5);
% 
% figure;hold
% set(gca,'fontsize',18)
% for ii = 1:length(CHIABV)
%     COL = (ii-1)/(length(CHIABV)-1)
%     CHI = CHIABV(ii);
% 
%     [CHIABS, KS, ~, ~] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
%     [EIG,EIGV] = gamma2_solvent(N,NM,LAM,FA,KV,[0, CHI; CHI, 0]*CHIABS,PHIP);
% 
%     plot(KV*RM, 1./EIG(:,1)/NM, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
%     plot(KV*RM, 1./EIG(:,2)/NM, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
% end
% 
% box on
% xlabel('Rq');
% ylabel('$\langle {\psi}^{\bf{\top}}(\vec{q}) {\psi}(-\vec{q}) \rangle/N_{m}v$',...
%     'Interpreter', 'Latex')
% xlim([min(KV), max(KV)]*RM)
% set(gca,'fontsize',18)
% set(gca,'xscale','log')
% set(gca,'yscale','log')