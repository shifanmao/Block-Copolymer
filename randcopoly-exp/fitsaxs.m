%% 
% This script uses analytical theory to fit experimental SAXS data of Polyimide-PEG 
% random copolymers.

clear;
MW = 1500;
WT = 30;
FILENAME = sprintf('PEG%dMW%d',WT,MW);

global sf qf N rm FA NFIT

% Load experiment data, SAXS data with q in A^(-1)
addpath('functions/')
data = load(strcat('exp-data/',FILENAME,'.csv'));  

% preset parameters in theory
N=100;  % total of 100 monomers
[FA,rm]=calcmol(MW,WT/100);

if (MW==1500)
    TV = [22,40:20:180];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
elseif (MW==900)
    TV = [22,40:20:160];  % temperature in degree C
    TK = TV+273.15;       % temperature in Kelvin
end

q = data(:,1);
s = data(:,2:end);
%% 
% Define a small q cutoff and estimate peak intensity S(q~=0).

if (MW==1500)
    qmin0 = 0.03803;
    qmin0f = 0.04094;
    qfinal = 0.22138;
elseif (MW==900)
    qmin0 = 0.05017;
    qmin0f = 0.08312;
    qfinal = 0.31367;
end

iqmin0 = find(q==qmin0);
iqmin0f = find(q==qmin0f);
iqfinal = find(q==qfinal);
smin = mean(s(iqmin0:iqmin0f,:));
sminstd = std(1./s(iqmin0:iqmin0f,:));

%% 
% Start fitting analytical S(q) to experimental data.
% number of fitting temperatures
NFIT = length(TV);
%NFIT = 6;

% define fitting range
qf = q(iqmin0:iqfinal);
sf = zeros(length(qf),NFIT);
for IT = 1:NFIT
	   sf(:,IT) = s(iqmin0:iqfinal,IT);
end

%% Fit with random copolymer model
% initial fit guess
% x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT

if strcmp(FILENAME,'PEG30MW1500')
    %x0 = [100,0.0,1.0,ones(1,NFIT)*2];
    %x0 = [114.0821    0.5147    0.1228    7.2675    7.5036    7.9280    8.1986    8.4047    8.5794    8.7466    8.8776    8.9516];
    %x0 = [190.4329    0.0040    1.0101    0.0112    0.0124    0.4639    0.9128    1.2622    1.5611 1.8442    2.0651    2.1917];
    %x0 = [191.4407   -0.0417    0.8458    0.0035    0.0167    0.8417    1.3773    1.7857    2.1328    2.4636    2.7224    2.8689];
    %x0 = [187.8140   -0.0310    0.6663    0.0022    0.0375    1.0726    1.7368    2.2428    2.6726    3.0824    3.4030    3.5845];
    %x0 = [184.5047   -0.0234    0.5201    0.0007    0.0775    1.3842    2.2188    2.8543    3.3939    3.9083    4.3108    4.5385];
    %x0 = [178.1441   -0.0107    0.2676    0.0001    0.3477    2.8104    4.3847    5.5836    6.6012    7.5723    8.3319    8.7609];
    %x0 = [183.0337   -0.1123    0.0350    0.0001   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000];
    %x0 = [181.8999   -0.0742    0.1698    0.0001    1.4735    5.4923    8.0641   10.0000   10.0000   10.0000   10.0000   10.0000];
    %x0 = [175.4019   -0.0132    0.2309    0.0001    0.7147    3.5288    5.3275    7.0984    7.9621    8.9557   10.0000   10.0000];
    x0 = [175.4019   -0.0132    0.2309    0.0001    0.7147    3.5288    5.3275    7.0984    7.9621    8.9557   10.0000   10.0000];
elseif strcmp(FILENAME,'PEG30MW900')
    %x0 = [500, 0.45,   0.08,  ones(1,NFIT)*10];
    x0 = [278.3541    0.3696    0.0930   10.1826   10.1926   10.4964   10.6912];
    %x0 = [273.2604    0.4732    0.1065    6.1926    6.2053    6.4995    6.6911    6.8765    7.1276    7.1276    7.2773];
elseif strcmp(FILENAME,'PEG35MW900')
    x0 = [260.5393    0.4182    0.0846    8.2573    8.3828    8.7970    8.9931    9.1637    9.2139    9.2317    9.3414];
else
	x0 = [500, 0.45,   0.08,  ones(1,NFIT)*10];
end

x0 = x0(1:NFIT+3);
lb = [10,   -0.2,  0.001, ones(1,NFIT)*1e-4];
ub = [1e5,   1.00,  10.0, ones(1,NFIT)*10.0];

% start fitting
options = optimset('MaxFunEvals',5,'Display','iter');
x = lsqnonlin(@saxsfitfunc,x0,lb,ub,options);
% x = x0;


%% Plot structure factors
figure;hold;set(gca,'fontsize',18)
for IT = 1:length(TV)
    if length(TV)==1
        col = 1;
    else
        col = (IT-1)/(length(TV)-1);
    end
    plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
end

% reconstruct fitted function
sfit = zeros(length(qf),NFIT);
for IT = 1:NFIT
    if NFIT==1
        col = 1;
    else
        col = (IT-1)/(length(TV)-1);
    end
    
    % mdl = @(x,qf) x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm)*x(3))';
    NM = x(3);
    R2=-0.5+0.5*exp(-2*NM)+NM;
    sfit(:,IT) = x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm/sqrt(R2)));
    sfit(:,IT) = sfit(:,IT)/NM;
    plot(qf*rm,sfit(:,IT),'--','LineWidth',2,'color',[col 0 1-col]) % theoretical fit
end

xlabel('R_Mq');ylabel('S(q)');
set(gca,'xscale','log');set(gca,'yscale','log');
%xlim([0.02,0.6]*rm);ylim([7,1e2]);
box on

% Save fitted results
SCALE = x(1);
LAM = x(2);
NM = x(3);
CHI = x(4:end);
save(strcat('savedata/',FILENAME,'.mat'),'SCALE','LAM','NM','CHI','TK','qf');

%% % Fit with Ornstein-Zernike function (Lorentzian)
% % initial fit guess
% % x(1) = intensity scale, x(2) = LAM, x(3) = NM, x(4:NFIT+4) = CHI_IT
% x0 = [1,ones(1,NFIT)*10];
% 
% % start fitting
% % options = optimset('TolX',1e-2,'TolFun',1e0,'MaxFunEvals',1e2,'MaxIter',1e2,...
% %     'Display','iter');
% % fit = lsqnonlin(@ozfit,x0,lb,ub,options);
% fit = lsqnonlin(@ozfit,x0,[],[]);
% 
% %% % reconstruct fitted function
% sfit = zeros(length(qf),NFIT);
% for IT = 1:NFIT
%     % mdl = @(x,qf) x(1)./(-2*x(3+IT)+s2invwlc(N,x(3),FA,x(2),qf*rm)*x(3))';
%     sfit(:,IT) = fit(1)*fit(1+IT).^2./(1+qf.^2*fit(1+IT)^2);
% end
% 
% figure;hold
% for IT = 1:NFIT
	     %     col = (IT-1)/(NFIT-1);
%     plot(q*rm,s(:,IT),'-','LineWidth',1,'color',[col 0 1-col]) % experiment data
%     plot(qf*rm,sfit(:,IT),'--','LineWidth',2,'color',[col 0 1-col]) % theoretical fit
% end
% xlabel('R_Mq');ylabel('S(q)');
% set(gca,'xscale','log');set(gca,'yscale','log');
% xlim([0.02,0.6]*rm);ylim([7,1e2]);
% box on
