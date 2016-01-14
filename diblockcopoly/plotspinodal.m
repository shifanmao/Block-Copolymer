% Calculate spinodal and critical wavelength
clear all;

FAV = 0.5;
N = 1e5;
NMV = logspace(-1,3,11);

% results to save
ks = zeros(length(NMV),1);
chis = zeros(length(NMV),1);

% start calcultion
inm = 1;
for NM = NMV
    [chis(inm),ks(inm)]=spinodal(N,NM,FAV);
    inm = inm+1;
end

figure;semilogx(NMV,chis*N)

R2=sqrt(r2(NMV))';
figure;semilogx(NMV,1./(ks.*R2))