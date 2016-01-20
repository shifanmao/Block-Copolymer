%% Calculate spinodal and critical wavelength
clear all;

FAV = 0.5;
NV = logspace(-2,3,21);

% results to save
ks = zeros(length(NV),1);
chis = zeros(length(NV),1);
d2gam2 = zeros(length(NV),1);

% start calcultion
inm = 1;
for N = NV
    [chis(inm),ks(inm),d2gam2(inm)]=spinodal(N,FAV);
    inm = inm+1;
end

figure;semilogx(NV,chis.*NV')
figure;loglog(NV,2*pi./ks)
figure;semilogx(NV,1./(ks.*sqrt(r2(NV))'))