function [x,sfit]=saxsfit(qf,sf,x0,N,NM,FA,LAM,CHIS)
% define fitting function
% param(1) = CHI
% param(2) = scale
% param(3) = rm
Sfun = @(param,qfit) param(2)./(-2*param(1)*CHIS+s2invwlc(N,NM,FA,LAM,qfit*param(3))*NM);

% start fitting
lb = [0.001,1,  1];
ub = [1.0,  1e3,1e3];

options = optimset('Display','on',...
    'TolX',1e-3,'TolFun',1e0,'MaxFunEvals',5e1,'MaxIter',1e1);
x = lsqcurvefit(Sfun,x0,qf,sf,lb,ub,options);

% evaluate fitted structure factor intensity
sfit = Sfun(x,qf);