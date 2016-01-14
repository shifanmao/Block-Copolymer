function [chis,ks,d2gam2]=spinodal(N,NM,FAV)
% spinodal.m :: This function predicts diblock copolymer phase transition spinodal
% and critical wavelength of quadratic density fluctuation in the well-mixed state
% Usage: [chis,ks,d2gam2]=spinodal(N,FAV)
% Parameters:
%   CHAIN, type of polymer (CHAIN=1, Guassian; =2, WLC; =3 Rigid Rod)
%   NM, number of Kuhn steps of total chain
%   FAV, range of A-type monomer fractions
%   ORDEig, number of eigenvalues
%   NumLayer, Number of residual layers
% Return:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation

% add paths
addpath('misc')
addpath('chainstats')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

% results to return
chis=zeros(length(FAV),1);     % spinodal
ks=zeros(length(FAV),1);        % critical wavelength of density fluctuations
d2gam2=zeros(length(FAV),1);   % inverse of susceptibility

for ii=1:length(FAV)
    FA=FAV(ii);
    fprintf('Step 1: Calculating spinodal at FA=%.2f, NM=%.2f\n',FA,NM)

    % find kstar
    G=@(k) gamma2(N,NM,FA,k,0);

    % initial guesses of kstar
    R2=r2(NM);
    k0=-1e-3;
    kf=1e3;
    ks(ii)=fminbnd(G,k0,kf);

    % spinodal
    chis(ii) = 0.5*G(ks(ii));

    % find susceptibility (second der. of gamma2 w/ k at kstar)
    dks = 1/sqrt(R2)*5e-2;

    if ks(ii)>1e-1  % central differences
        d2gam2(ii) = (G(ks(ii)+dks)-2*G(ks(ii))+G(ks(ii)-dks))/(dks^2);
    else  % forward differences
        d2gam2(ii) = (G(ks(ii)+2*dks)-2*G(ks(ii)+dks)+G(ks(ii)))/(dks^2);
    end
    d2gam2(ii) = -1/(N*R2)*d2gam2(ii)./power(G(ks(ii)),2);
end