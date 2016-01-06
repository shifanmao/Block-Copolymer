function [chis,ks,d2gam2]=spinodal(M,NM,LAM,FAV)
%% spinodal.m :: This function predicts phase transition spinodal
% Usage: [chis,ks,d2gam2]=spinodal(M,NM,LAM,FAV)
% and critical wavelength of quadratic density fluctuation in the well-mixed state
% Parameters:
%   NM, number of Kuhn steps of total chain
%   FAV, range of A-type monomer fractions
%   ORDEig, number of eigenvalues
%   NumLayer, Number of residual layers
% Return:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation
%   d2gam2, second derivative of density-density correlation

% import directory
addpath('chainstats')

% results to return
chis=zeros(length(FAV),1);     % spinodal
ks=zeros(length(FAV),1);       % critical wavelength of density fluctuations
d2gam2=zeros(length(FAV),1);   % inverse of susceptibility

for ii=1:length(FAV)
FA=FAV(ii);
disp(['Step 1: Calculating spinodal at FA=',num2str(FA), ', NM=',num2str(NM),', LAM=',num2str(LAM)])

%% find kstar
G=@(k) gamma2(M,NM,LAM,FA,k,0);

% initial guesses of kstar
R2 = r2(NM);
k0=-1e-2/(sqrt(R2));
kf=1e1/(sqrt(R2));

options = optimset('Display','off',...
    'TolX',1e-5,'TolFun',1e-5,'MaxFunEvals',1e5,'MaxIter',1e5);
ks(ii) = fminbnd(G,k0,kf,options);
if ks(ii)*sqrt(R2)<=1e-1
    ks(ii)=0;
end

%% find spinodal
chis(ii) = 0.5*G(ks(ii));

%% find susceptibility (second der. of gamma2 w/ k at kstar)
dks = 1/sqrt(R2)*5e-2;

if ks(ii)>1e-1  % central differences
    d2gam2(ii) = (G(ks(ii)+dks)-2*G(ks(ii))+G(ks(ii)-dks))/(dks^2);
else  % forward differences
    d2gam2(ii) = (G(ks(ii)+2*dks)-2*G(ks(ii)+dks)+G(ks(ii)))/(dks^2);
end
d2gam2(ii) = -1/(NM*R2)*d2gam2(ii)./power(G(ks(ii)),2);
end