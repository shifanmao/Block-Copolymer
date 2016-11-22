# Chainstats
This folder contains functions that calculate structure factors of random copolymers (2,3,and 4-point correlations).

Example 1: Two-point correlations of Gaussian random copolymer
`
clear;

N = 100;  % Number of Monomers
NM = 100; % Number of Kuhn steps per monomer
LAM = 0.0;  % Degree of chemical correlation
FA = 0.5; % Fraction of A-type monomers
kv = logspace(-4,2,100);  % Wavevector

s2 = zeros(length(kv),2,2);
SAA = zeros(length(kv),1);
SAB = zeros(length(kv),1);
SBA = zeros(length(kv),1);
SBB = zeros(length(kv),1);
for ii = 1:length(kv)
  k = kv(ii);  % wavevector
  s2(ii,:,:)=s2gc(N,NM,LAM,FA,k)/power(NM*N,2);
  SAA(ii) = s2(ii,1,1);
  SAB(ii) = s2(ii,1,2);
  SBA(ii) = s2(ii,2,1);
  SBB(ii) = s2(ii,2,2);
end

figure;
semilogx(kv,s2(:,1,1));  % Plots SAA for different k

% check normaliztion
figure;
semilogx(kv,SAA+SAB+SBA+SBB);  % Plots SAA+SAB+SBA+SBB for different k (expect to go to 1 at zero k)

% plot melt structure factor
CHI = 0;
S = SAA+SAB+SBA+SBB;
W = SAA.*SBB-SAB.*SBA;
Gamma2 = -2*CHI+S./W;
figure;semilogx(kv,1./Gamma2)

% calculate solvent structure factor
...
`