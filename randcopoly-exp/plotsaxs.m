%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;
MW = 900;
WT = 40;
FILENAME = sprintf('PEG%dMW%d',WT,MW);

% Load experiment data.
addpath('functions/')
%load(strcat('savedata/',FILENAME,'.mat'));
data = load(strcat('exp-data/',FILENAME,'.csv'));

% preset parameters in theory
N=100;  % total of 100 monomers
[FA,rm]=calcmol(1500,0.3);

q = data(:,1);
s = data(:,2:end);

% extract domain with low background noise
% iq0 = find(q==qf(1));iqf = find(q==qf(end));
% qf = q(iq0:iqf);
% sf = s(iq0:iqf,:);

if (MW==1500)
    TV = [22,40:20:180];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
    qmin0 = 0.03803;
elseif (MW==900)
    TV = [22,40:20:160];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
    qmin0 = 0.05017;
end

%% PLOT 1: Structure factors
NFIT = length(TK);

figure;hold
set(gca,'fontsize',18)
for IT = 1:NFIT
   col = (IT-1)/(NFIT-1);
   plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
end

% for IT = 1:NFIT
%    col = (IT-1)/(NFIT-1);
%    plot(qf*rm,sf(:,IT),'+-','LineWidth',1,'color',[col 0 1-col]) % experiment data
% end

% plot process
xlabel('R_Mq');ylabel('S(q)');
set(gca,'xscale','log');set(gca,'yscale','log');
ylim([1e1,1e3]);
plot([1.35,1.35],[10,1e3],'k--')
box on

% save to data
% SAVEFILENAME = strcat('exp-data/','SEXP_',FILENAME);
% dlmwrite(SAVEFILENAME,[qf*rm,sf]);