%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;
%close all
MW = 1500;
WT = 30;
FILENAME = sprintf('PEG%dMW%d',WT,MW);

% Load experiment data.
addpath('functions/')
data = load(strcat('exp-data/',FILENAME,'.csv'));

% preset parameters in theory
N=100;  % total of 100 monomers
[FA,rm]=calcmol(MW,WT/100);
q = data(:,1);
s = data(:,2:end);
load(strcat('savedata/',FILENAME,'.mat'));

%load(strcat('savedata/PEG35MW900.mat'));
if strcmp(FILENAME,'PEG30MW900')
    TK = TK(1:end-1);
    qmin0 = 0.0577;
    iq = find(q==qmin0);
    SCALE = SCALE*0.;
else
    iq = find(q==qf(1));
end
iq = find(q==qf(1));

%x0 = [114.0821    0.5147    0.1228  7.2675    7.5036    7.9280    8.1986    8.4047    8.5794    8.7466    8.8776    8.9516];
%x0 = [190.4329    0.0040    1.0101    0.0112    0.0124    0.4639    0.9128    1.2622    1.5611 1.8442    2.0651    2.1917];
%x0 = [191.4407   -0.0417    0.8458    0.0035    0.0167    0.8417    1.3773    1.7857    2.1328    2.4636    2.7224    2.8689];
%x0 = [183.0337   -0.1123    0.0350    0.0001   10.0000];
%x0 = [177.2075   -0.0090    0.2344    0.0001    0.4391    3.2386    5.0280    6.3911    7.5480    8.6527    9.5175   10.0000];
x0 = [175.4019   -0.0132    0.2309    0.0001    0.7147    3.5288    5.3275    7.0984    7.9621    8.9557   10.0000   10.0000];

x = x0;
% Save fitted results
SCALE = x(1);
LAM = x(2);
NM = x(3);
CHI = x(4:end);
R2=-0.5+0.5*exp(-2*NM)+NM;
sqrt(R2)

%% PLOT 1: Structure factors
NFIT = length(TK);
sfit = zeros(length(qf),NFIT);
for IT = 1:NFIT
   sfit(:,IT) = SCALE./(-2*CHI(IT)+s2invwlc(N,NM,FA,LAM,qf*rm/sqrt(R2)));
   %sfit(:,IT) = SCALE./(-2*CHI(IT)+s2invwlc(N,NM,FA,LAM,qf*rm));
end

figure;hold
set(gca,'fontsize',18)
for IT = 1:NFIT
   col = (IT-1)/(NFIT-1);
   plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
   plot(qf*rm,sfit(:,IT)/NM,'--','LineWidth',2,'color',[col 0 1-col]) % theoretical fit
end
xlabel('R_Mq');ylabel('S(q)');
set(gca,'xscale','log');set(gca,'yscale','log');
% xlim([0.02,0.6]*rm);ylim([7,1e2]);
box on

%% PLOT 2: CHI vs T, and 1/S(q*) vs T
figure;
subplot(2,1,1);hold;set(gca,'fontsize',18)
for IT = 1:NFIT
    col = (IT-1)/(NFIT-1);
    plot(-1./TK(IT),1./s(iq,IT),'.','color',[col 0 1-col],...
        'MarkerSize',25);
end
[xp,yp] = plotregressor(-1./TK,1./s(iq,:),1:3,1:length(TK));
plot(xp,yp,'k--','LineWidth',2);
xlabel('-1/T (K^{-1})');ylabel('S^{-1}(q=0)');box on

subplot(2,1,2);hold;set(gca,'fontsize',18)
for IT = 1:NFIT
    col = (IT-1)/(NFIT-1);
    plot(-1./TK(IT),CHI(IT)*NM,'.','color',[col 0 1-col],...
        'MarkerSize',25);
end
[xp,yp] = plotregressor(-1./TK,CHI*NM,1:NFIT,1:NFIT);
plot(xp,yp,'k--','LineWidth',2);
%ylim([0.8,1.4])
xlabel('-1/T (K^{-1})');ylabel('\chiN_M');box on

% %% PLOT 3: CHI vs 1/S(q*)
% figure;hold;set(gca,'fontsize',20)
% plot(CHI(:)*NM,NM./s(iq,:),'k-','linewidth',1.2);
% plot(CHI(:)*NM,NM./sfit(1,:),'k-','linewidth',1.2);
% for IT = 1:NFIT
%     col = (IT-1)/(NFIT-1);
%     plot(CHI(IT)*NM,NM./s(iq,IT),'o-','color',[col 0 1-col],...
%         'MarkerSize',10,'linewidth',2);
% %     plot(CHI(IT)*NM,NM./sfit(1,IT),'x-','color',[col 0 1-col],...
% %         'MarkerSize',10,'linewidth',2);
% end
% xlabel('\chiN_M');ylabel('S^{-1}(q=0)');box on

%% PLOT 4: Phase diagram
FAV = linspace(0.1,0.9,51);
CHISV = zeros(length(FAV),1);
for ii = 1:length(FAV)
   [kval,sval,d2gam2]=kmaxwlc(N,NM,FAV(ii),LAM);
   CHISV(ii)=0.5*sval;  % spinodal
end

figure;hold;set(gca,'fontsize',20)
plot(FAV,CHISV*NM,'k-','linewidth',2)
for IT=1:NFIT
	col = (IT-1)/(NFIT-1);
    plot(FA,CHI(IT)*NM,'.','color',[col 0 1-col],...
         'LineWidth',2,'MarkerSize',25);
    plot(0.2382,CHI(IT)*NM,'.','color',[col 0 1-col],...
         'LineWidth',2,'MarkerSize',25);
end
xlabel('PEG mol%');ylabel('\chi_S v N_M')
%ylim([0.5,1.5]);xlim([0.1,0.9]);box on