% test the integral considered by Stepanow
clear;
%close all
addpath('../');
addpath('../eigcalc/');
addpath('../integrals/');
addpath('../../misc/');
addpath('../../functions/');

N = 1e0;
FA = 0.5;
[chis,ks,d2gam2]=spinodal(N,0.5);

RM = sqrt(N/6);
QM = linspace(1e-2,10,200)'/RM;
gam4 = zeros(length(QM),1);
gam2 = zeros(length(QM),1);
QS=[1,0,0]'*ks;
for ii = 1:length(QM)
    ii
    Q1=[1,0,0]'*QM(ii);
    gam4(ii) = gamma4allq(N,FA,Q1,-Q1,QS,-QS);
end

% CHI = 10.3/N;
CHI = 6.0/N;
for ii = 1:length(QM)
    gam2(ii) = gamma2(N,FA,QM(ii),CHI);
end
gam4star = gamma4allq(N,FA,QS,-QS,QS,-QS);

int = (QM*RM).^2.*gam4./gam2;
braz = (ks*RM).^2*gam4star./gam2;

figure;hold;set(gca,'fontsize',20)
plot(QM*RM,int,'k-','linewidth',2)
plot(QM*RM,braz,'k--','linewidth',2)
xlabel('qRg');ylabel('Integrand');box on