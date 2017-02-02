clear
close all
addpath(genpath('chainstats/'))
addpath('figures/')

%%% PARAMETERS %%%
N = 100;
NM = 0.1;
FA = 0.3;
PHIP = 0.5;   % Fraction of polymers

% WAVEVECTOR
NK = 100;
RM = (r2(NM))^0.5; % Normalization factor
KV = logspace(-2,2,NK)/RM;  % Wavevector

% CHI PARAMETERS
CHIBS = 0; % Flory-Huggins factor between B and S
CHIAS = 0;

%% PLOT 1: Structual factor plot
CHIABV = [0:0.2:0.8];
for LAM = [-0.75,0.00]
    fprintf(' Calculating Structure Factor at lambda=%.2f\n', LAM)
    
    figure;hold;set(gca,'fontsize',18)
    CHIABS = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
    for ii = 1:length(CHIABV)
        COL = (ii-1)/(length(CHIABV)-1);
        CHIAB = CHIABV(ii)*CHIABS; % Flory-Huggins factor between A and S
        CHIBA = CHIAB;
        CHI = [CHIAS,CHIAB;CHIBA,CHIBS];

        [EIG1,EIG2,~,~,~,~]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,PHIP);
        loglog(KV*RM, 1./EIG1/NM, '-', 'Linewidth', 1.5,'color',[COL 0 1-COL]);
        loglog(KV*RM, 1./EIG2/NM, '--', 'Linewidth', 1.5,'color',[COL 0 1-COL]);
    end
    xlabel('R_M\cdotq')
    ylabel('$\left<\psi(q)\psi(-q)\right>$', 'interpreter', 'latex')
    title(strcat('$\phi_{P}$', sprintf(' = %.2f', PHIP), ...
                 '$, f_{A}$', sprintf(' = %.2f', FA), ...
                 ', $\lambda$', sprintf(' = %.2f', LAM)), ...
                 'interpreter', 'latex')
    set(gca,'XScale','log', 'YScale', 'log')
    box on
    ylim([1e-3,1e1])
end

% %% PLOT 2: Critical wavemode direction
% for LAM = [-0.75,0.00]
%     fprintf(' Calculating Critical Wavevector at lambda=%.2f\n', LAM)
%     
%     f = figure('position', [0, 0, 600, 550]);hold
%     t = linspace(0,2*pi,100);x = sin(t);y = cos(t);
%     set(gca,'fontsize',18)
%     plot(x,y,'k-','linewidth',2)
%     title(strcat('$\phi_{P}$', sprintf(' = %.2f', PHIP), ...
%                  '$, f_{A}$', sprintf(' = %.2f', FA), ...
%                  ', $\lambda$', sprintf(' = %.2f', LAM)), ...
%                  'interpreter', 'latex')
%     box on
% 
%     CHIABS = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
%     CHIABV = [0:0.2:0.8]*CHIABS;
%     for ii = 1:length(CHIABV)
%         COL = (ii-1)/(length(CHIABV)-1);
%         CHIAB = CHIABV(ii); % Flory-Huggins factor between A and S
%         CHIBA = CHIAB;
%         CHI = [CHIAS,CHIAB;CHIBA,CHIBS];
%         
%         [~,~,EIGV1,EIGV2,KS1,KS2]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,PHIP);
%         
%         KS1IND = find(KV>=KS1,1);
%         KS2IND = find(KV>=KS2,1);
%         plot([0, -EIGV1(KS1IND,1)], [0, -EIGV1(KS1IND,2)],'-','linewidth',2,...
%               'color',[COL 0 1-COL]);
%         plot(-EIGV1(KS1IND,1), -EIGV1(KS1IND,2),'.','markersize',20,'color',[COL 0 1-COL]);
%         plot([0, -EIGV2(KS2IND,1)], [0, -EIGV2(KS2IND,2)],'--','linewidth',2,...
%               'color',[COL 0 1-COL]);
%         plot(-EIGV2(KS2IND,1), -EIGV2(KS2IND,2),'.','markersize',20,'color',[COL 0 1-COL]);
%     end
% end
% 
% %% PLOT 3: Phase diagram: chis*NM vs LAM
% LAMV = linspace(-1,1,100);
% CHIS = zeros(length(LAMV),1);
% for ii = 1:length(LAMV)
%     fprintf(' Calculating Spinodal at lambda=%.2f\n', LAM)
%     LAM = LAMV(ii);
%     CHIS(ii) = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
% end
% 
% figure;set(gca,'fontsize',18)
% plot(LAMV,CHIS*NM*4*FA*(1-FA),'k-','linewidth',2)
% xlabel('$\lambda$', 'interpreter', 'latex')
% ylabel('$4f_{A}f_{B}v\chi_{\mathrm{S}}N_{\mathrm{M}}$', 'interpreter', 'latex')

% %% PLOT 4: Critical wavemode k* vs. LAM and CHI% LAMV = linspace(-1,1,20);
% 
% LAMV = linspace(-1,1,100);
% CHIABV = [0:0.2:0.8];
% 
% figure;hold;set(gca,'fontsize',18)
% for ii = 1:length(CHIABV)
%     fprintf(' Calculating Spinodal at CHI=%.2fCHIS\n', CHIABV(ii))
% 
%     COL = (ii-1)/(length(CHIABV)-1);
%     CHIABS = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
%     CHIAB = CHIABV(ii)*CHIABS; % Flory-Huggins factor between A and S
%     CHIBA = CHIAB;
%     CHI = [CHIAS,CHIAB;CHIBA,CHIBS];
% 
%     KS1 = zeros(length(LAMV),1);
%     KS2 = zeros(length(LAMV),1);
%     NK = 1000;
%     KV = linspace(0.1,5,NK)/RM;  % Wavevector
%     for jj = 1:length(LAMV)
%         LAM = LAMV(jj);
%         [~,~,~,~,KS1(jj),KS2(jj)]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,PHIP);
%     end
%     KS1(KS1*RM<=0.11) = 0;
%     KS2(KS2*RM<=0.11) = 0;
%     
%     plot(LAMV,KS1*RM,'-','linewidth',2,'color',[COL 0 1-COL]);
%     plot(LAMV,KS2*RM,'--','linewidth',2,'color',[COL 0 1-COL]);
% end
% box on
% xlabel('$\lambda$', 'interpreter', 'latex')
% ylabel('$R_{M}q^{*}$', 'interpreter', 'latex')

%% Save figures
filename = sprintf('FA%.2fPHIP%.2fNM%.2f',FA,PHIP,NM);
figure(1);saveas(gca,strcat('figures/S_',filename,'lamnp75.eps'),'epsc')
figure(2);saveas(gca,strcat('figures/S_',filename,'lam0.eps'),'epsc')
% figure(3);exportfig(strcat('figures/psi_',filename,'lamnp75.jpg'))
% figure(4);exportfig(strcat('figures/psi_',filename,'lam0.jpg'))
% figure(5);saveas(gca,strcat('figures/chis_',filename,'.eps'),'epsc')
% figure(6);saveas(gca,strcat('figures/ks_',filename,'.eps'),'epsc')