function [KS,SINV,D2S,ERR]=calcserr(K,S,G)
% use Lorentzian fit to estimate peak location, intensity and errors

IND = find(S==max(S));IND = IND(1);
KS = K(IND);
SINV = 1/S(IND);

NUMFIT = 4;
if IND>NUMFIT  % central differences
    Kfit = K(IND-NUMFIT:IND+NUMFIT);
    Sfit = S(IND-NUMFIT:IND+NUMFIT);

    % local fit to Lorentzian (three parameter fit)
    options = optimset('Display','off',...
       'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e10,'MaxIter',1e10);
    x0 = [KS,SINV,1];   % x = [ks_sim,sinv_sim,d2s]
    lb = [0.1,0.01,0.1]; ub = [5,10,10];
    fun = @(x,Kfit) x(2) + (1/2)*(1/G)*x(3)*(Kfit-x(1)).^2;
    [xhat,~,resid,~,~,~,J] = lsqcurvefit(fun,x0,Kfit,1./Sfit,lb,ub,options);
    [~,sehat] = nlparcinew(xhat,resid,'jacobian',J);      % %estimate error from the fit
    x=xhat;se=sehat;
elseif IND<=NUMFIT  % forward differences
    Kfit = K(IND:IND+NUMFIT);
    Sfit = S(IND:IND+NUMFIT);

    % local fit to Lorentzian with fixed peak at zero k* (two parameter fit)
    options = optimset('Display','off',...
       'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e10,'MaxIter',1e10);
    x0 = [SINV,1];   % x = [sinv_sim,d2s]
    lb = [0.01,0.1]; ub = [10,10];
    fun = @(x,Kfit) x(1) + (1/2)*(1/G)*x(2)*(Kfit).^2;
    [xhat,~,resid,~,~,~,J] = lsqcurvefit(fun,x0,Kfit,1./Sfit,lb,ub,options);
    [~,sehat] = nlparcinew(xhat,resid,'jacobian',J);      % %estimate error from the fit
    x=[0,xhat];se=[0,sehat'];
else
    error('peak location has large variability')
end

KS = x(1);
SINV = x(2);
D2S = x(3);
ERR = se(1:3);
