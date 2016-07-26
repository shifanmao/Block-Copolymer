%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;close all
FILENAME = 'PEG40MW1500';
paramcp = 6; chemcp = 1;
FOLDERNAME = 'scalcbatch-05-31-16/sdata-6-1/';

%% Load experiment data.
addpath('functions/')
data = load(strcat('exp-data/',FILENAME,'.csv'));  % SAXS data with q in A^(-1)
                                           % 30wt ~= 16mol%
% preset parameters in theory
N=100;  % total of 100 monomers
q = data(:,1);
s = data(:,2:end);

if ( strcmp(FILENAME,'PEG30MW900') || strcmp(FILENAME,'PEG30MW900'))
    load(strcat('savedata/PEG30MW1500.mat'));
    TV = [22,40:20:160];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
    [FA,rm]=calcmol(900,0.3);
elseif ( strcmp(FILENAME,'PEG30MW1500') || strcmp(FILENAME,'PEG40MW1500') )
    load(strcat('savedata/PEG30MW1500.mat'));
    TV = [22,40:20:160];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
    [FA,rm]=calcmol(1500,0.3);
end

%% PLOT 1: Structure factors
figure;hold;
set(gca,'fontsize',18)
NFIT = length(TK);
[kval,sval,d2gam2]=kmaxwlc(N,NM,0.2382,LAM);
CHIS=0.5*sval;  % spinodal
s = fliplr(s);CHI = fliplr(CHI);
leg = {};
for IT = 1:NFIT
    col = 1-(IT-1)/(NFIT-1);
    plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
    leg = [leg,strcat('Exp. \chi/\chi_S',sprintf(' = %.2f',CHI(IT)/CHIS))];
end

%% Load Simulation data.
SCALEMC = 27;
CHIMC = load(strcat(FOLDERNAME,'Sdata/chilist'));CHIMCS = CHIMC(11);
CHIMC = flipud(CHIMC([20,22]));
SMC = zeros(41,2,length(CHIMC));
for ICHI = 1:length(CHIMC)
    SMC(:,:,ICHI) = load(strcat(FOLDERNAME,...
        sprintf('Sdata/SMC_SIM%dCHEM%dCHI%.8f',paramcp,chemcp,CHIMC(ICHI))));
end
for ICHI = 1:length(CHIMC)
    col = (ICHI-1)/(length(CHIMC)-1);
    plot(SMC(:,1,ICHI),SMC(:,2,ICHI)*SCALEMC,'s-',...
        'LineWidth',2,'color',[col 0 1-col])
    leg = [leg,strcat('Sim. \chi/\chi_S',sprintf(' = %.2f',CHIMC(ICHI)/CHIMCS))];
end


CHIMC = load(strcat(FOLDERNAME,'Sdata/chilist'));CHIMCS = CHIMC(11);
CHIMC = flipud(CHIMC([9,10,11]));
SMC = zeros(41,2,length(CHIMC));
for ICHI = 1:length(CHIMC)
    SMC(:,:,ICHI) = load(strcat(FOLDERNAME,...
        sprintf('Sdata/SMC_SIM%dCHEM%dCHI%.8f',paramcp,chemcp,CHIMC(ICHI))));
end

for ICHI = 1:length(CHIMC)
    col = 1-(ICHI-1)/(length(CHIMC)-1);
    plot(SMC(:,1,ICHI),SMC(:,2,ICHI)*SCALEMC,'o--',...
        'LineWidth',2,'color',[col 0 1-col])
    leg = [leg,strcat('Sim. \chi/\chi_S',sprintf(' = %.2f',CHIMC(ICHI)/CHIMCS+0.2))];
end

set(gca,'xscale','log');set(gca,'yscale','log');
legend(leg)

xlabel('R_Mq');ylabel('S(q)');
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([0.7,11]);ylim([5,6e2]);
box on