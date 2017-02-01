clear
addpath(genpath('chainstats/'))

%%% PARAMETERS %%%
N = 100;
NM = 0.1;
LAM = 0;
FA = 0.5;

% VOLUME FRACTIONS
FP = 0.5;   % Fraction of polymers

NK = 100;
RM = (r2(NM))^0.5; % Normalization factor
KV = logspace(-2,2,NK)/RM;  % Wavevector

% CHI PARAMETERS
CHIBS = 0; % Flory-Huggins factor between B and S
CHIAS = 0;

%% Structual factor plot
figure;hold
CHIABS = solvnt_spin(N,NM,LAM,FA,KV,FP);
CHIABV = [0:0.2:0.8]*CHIABS;
for ii = 1:length(CHIABV)
    COL = (ii-1)/(length(CHIABV)-1);
    CHIAB = CHIABV(ii); % Flory-Huggins factor between A and S
    CHIBA = CHIAB;
    CHI = [CHIAS,CHIAB;CHIBA,CHIBS];
    
    [EIG1,EIG2,EIGV1,EIGV2,THETA1,THETA2,KS1,KS2]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,FP);
    loglog(KV*RM, 1./EIG1, '-', 'Linewidth', 1.5,'color',[COL 0 1-COL]);
    loglog(KV*RM, 1./EIG2, '--', 'Linewidth', 1.5,'color',[COL 0 1-COL]);
end
xlabel('R_M\cdotq')
ylabel('<\psi_A(q)\psi_B(-q)>')
title(sprintf('FP = %f', FP))
set(gca,'XScale','log', 'YScale', 'log')
box on
% 
% %% Phase diagram: chis*NM vs LAM
% LAMV = linspace(-0.99,0.99,20);
% CHIS = zeros(length(LAMV),1);
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii)
%     CHIS(ii) = solvnt_spin(N,NM,LAM,FA,KV,FP);
% end
% 
% figure;
% plot(LAMV,CHIS*NM)