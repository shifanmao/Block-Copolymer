function [chis,ks,d2gam2]=spinodal_solvnt(N,NMV,LAM,FAV,FP)
% spinodal.m :: This function predicts diblock copolymer phase transition spinodal
% and critical wavelength of quadratic density fluctuation in the well-mixed state
% Usage: [chis,ks,d2gam2]=spinodal(N,FAV)
% Inputs:
%   NMV, number of statistical steps of total chain
%   FAV, range of A-type monomer fractions
% Outputs:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation
%   d2gam2, second derivative of structure factor around peak

% results to return
chis=zeros(length(FAV),length(NMV));     % spinodal
ks=zeros(length(FAV),length(NMV));        % critical wavelength of density fluctuations
d2gam2=zeros(length(FAV),length(NMV));   % inverse of susceptibility

CHIBS = 0;
CHIAS = 0;
for ii=1:length(FAV)
    FA=FAV(ii);
    for jj=1:length(NMV)
        NM=NMV(jj);
        fprintf('Step 1: Calculating spinodal at N=%.2e,FA=%.2f\n',NM,FA)

        NK = 100;
        RM = (r2(NM))^0.5; % Normalization factor
        KV = logspace(-2,2,NK)/RM;  % Wavevector
        for CHIAB = linspace(0,10,100);
            CHIBA = CHIAB;
            CHI = [CHIAS,CHIAB;CHIBA,CHIBS];
            
            % find kstar
            [~,~,~,~,~,~,KS1,~]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,FP);
            G = gamma2_solvent(N,NM,LAM,FA,KS1,CHI,FP);
            
            if G<2e-2
                ks(ii,jj) = KS1;
                chis(ii,jj) = CHIAB;
                break
            end
        end
    end
end