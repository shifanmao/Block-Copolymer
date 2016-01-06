%% Calculate spinodal and critical wavelength
clear all;

FAV = 0.5;
NV = logspace(-1,3,21);

% results to save
ks = zeros(length(NV),1);
chis = zeros(length(NV),1);

% start calcultion
inm = 1;
for N = NV
    [chis(inm),ks(inm)]=spinodal(N,FAV);
    inm = inm+1;
end

figure;semilogx(NV,chis.*NV')