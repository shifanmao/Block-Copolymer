clear;

N = 100;
FA = .2;
NM = .01;
LAM = -0.22;

RM = sqrt(r2(NM));
KV = logspace(-2, 2, 100)/RM;

S2E = zeros(1, length(KV));
S2R = zeros(1, length(KV));

for ii = 1:length(KV);
    k = KV(ii);
    S2 = s2wlc_explicit(N,NM,LAM,FA,k);
    S2E(ii) = S2(1, 2);
    
    S2 = s2wlc_resum(N,NM,LAM,FA,k);
    S2R(ii) = S2(1, 2);
end

figure;
plot(KV*RM, S2E, '+', KV*RM, S2R,'-')
set(gca,'xscale','log');set(gca,'yscale','log')