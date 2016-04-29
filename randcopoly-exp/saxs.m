clear;
close all
addpath('functions/')

% load experiment data
TV = [22,40:20:180];  % temperature in degree C
TK = TV+273.15;       % temperature in Kelvin
data = load('exp-data/PEG30undoped.csv');  % SAXS data with q in A^(-1)
                                           % 30wt ~= 16mol%
q = data(:,1);
s = data(:,2:end);

% define small q cutoff
% % if PEG wt% = 0.30
% qmin0 = 0.03803;
% qminf = 0.04094;

% if PEG wt% = 0.40
qmin0 = 0.04155;
qminf = 0.04472;

iqmin0 = find(q==qmin0);
iqminf = find(q==qminf);
smin = mean(s(iqmin0:iqminf,:));
sminstd = std(1./s(iqmin0:iqminf,:));

% find temperature vs. peak intensity
figure;
hold;
xp = -1./TK;
yp = power(smin,-1);
errorbar(xp,yp,sminstd,'ko')
xlabel('-T^{-1} K^{-1})');ylabel('S^{-1}(q=0)')

% linear regression at q*=0
xfit = 1:3;  % linear fit regime
x = xp(xfit)';
X = [ones(length(x),1) x];
y = yp(xfit)';
b = X\y;
xp = linspace(-4,-1,100)*1e-3;
yp = b(2)*xp+b(1);
plot(xp,yp,'k--');
ylim([0,.1])

% analytical structure factor
N=100;  % total of 100 monomers
NM=0.1; % each monomer has 0.1 Kuhn steps
LAM=0.00; % ideal random copolymer
FA=0.165;    % equal chemical composition

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-2;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=101;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

%% find chi paramatrized with T
% CHI = A-B/T
TS = b(2)/b(1);
B = -b(2)/2;
A = CHIS*NM+B/TS;

I=9;
rm = 32.05; % end-to-end distance of a monomer in unit Angstrom
rm = 1;
figure;hold
plot(q*rm,s(:,I));






%% Make a plot
% at which temperature?
I = 9;

% normalization factor
CHI0 = 0.5;
Smin = power(-2*CHI0+1/(FA*(1-FA)),-1);  % zero q analytical structure factor

rm = 32.05; % end-to-end distance of a monomer in unit Angstrom
snorm = Smin/smin(I)*0.8;

% evaluate s2inv
[SINV]=s2invwlc(N,NM,FA,LAM,K);
[KS,SS]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*SS;

figure;hold
plot(q*rm,s(:,I));
plot([qmin0,qmin0]*rm,[7,1e2],'k-','linewidth',3)
plot([qminf,qminf]*rm,[7,1e2],'k-','linewidth',3)
loglog(RM*K,1./(-2*CHI0+SINV*NM),'-','LineWidth',2)

xlabel('R_Mq');ylabel('S(q)')
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([0.02,0.6]*rm);ylim([7,1e2])

%% find a region to fit
indmin = find(q*rm >= 2.5);indmin = indmin(1);
indmax = find(q*rm >= 21.5);indmax = indmax(1);
qmin = q(indmin);
qmax = q(indmax);
plot([qmin*rm,qmin*rm],[1e-3,1e3],'k--')
plot([qmax*rm,qmax*rm],[1e-3,1e3],'k--')
qfit = q(indmin:indmax);
sfit = s(indmin:indmax,I);

% define fitting function
% param(1) = CHI
% param(2) = scale
% param(3) = rm
Sfun = @(param,K) param(2)./(-2*param(1)*CHIS+s2invwlc(N,NM,FA,LAM,K*param(3))*NM);

% start fitting
% x0 = [0.01,5e2,500];
x0 = [0.001,snorm,rm];
lb = [0.001,1,  1];
ub = [1.0,  1e3,1e3];

options = optimset('Display','on',...
    'TolX',1e-3,'TolFun',1e0,'MaxFunEvals',1e1,'MaxIter',1e1);
% x = lsqcurvefit(Sfun,x0,qfit,sfit)
x = lsqcurvefit(Sfun,x0,qfit,sfit,lb,ub,options)
plot(q(indmin:indmax),Sfun(x,q(indmin:indmax)),'r--','LineWidth',2)

loglog(qfit,Sfun(x0,qfit),'k--')
% plot(RM*K,Sfun(x,K),'-','LineWidth',2)

xlabel('R_Mq');ylabel('S(q)')
set(gca,'xscale','log');set(gca,'yscale','log');


















% make a plot
figure;hold
for T = TV
    ii = find(TV==T);
    
    % current SAXS
    sp = s(:,ii);
    plot(q*rm,sp*snorm);
end
plot([qmin0,qmin0]*rm,[7,1e2]*snorm,'k-','linewidth',3)
plot([qminf,qminf]*rm,[7,1e2]*snorm,'k-','linewidth',3)
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([0.02,0.6]*rm);ylim([7,1e2]*snorm)

figure;hold;
I=1;
plot(q*rm,s(:,1)*snorm);
% plot(RM*K,1./(-2*CHI(I)+SINV))
set(gca,'xscale','log');set(gca,'yscale','log');
xlim([0.02,0.6]*rm);ylim([7,1e2]*snorm)

% analytical structure factor
addpath('functions/');
N=100;  % total of 100 monomers
NM=0.1; % each monomer has 0.1 Kuhn steps
LAM=0.00; % ideal random copolymer
FA=0.16;    % equal chemical composition

% find spinodal CHIS
[kval,sval]=kmaxwlc(N,NM,FA,LAM);
CHIS=0.5*sval;
CHI=CHIS*[0 0.2 0.4 0.6 0.8];

RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
K0=1e-2;  % minimum wavevector
KF=1e2;   % maximum wavevector
NK=201;  % number of wavevectors
K=transpose(logspace(log10(K0),log10(KF),NK))/RM;

% evaluate s2inv
[SINV]=s2invwlc(N,NM,FA,LAM,K);
figure;hold
for I=1:length(CHI)
    COL=(I-1)/(length(CHI)-1);
    loglog(RM*K,1./(-2*CHI(I)+SINV),'-','LineWidth',2,'Color',[COL 0 1-COL])
end
xlabel('R_Mq');ylabel('S(q)')
set(gca,'xscale','log');set(gca,'yscale','log');
axis([K0 KF 1e-4 1e-1]);box on