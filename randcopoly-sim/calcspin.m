%% Calculate spinodal and critical wavelength
clear all;

M=100;  % number of monomers
FA=0.5;  % fraction of A monomers
LAMV = linspace(-1,1-1e-3,51);  % degree of chemical correlation
% NMV = logspace(-2,2,11);  % number of Kuhn steps per monomer
NMV=1e6;

% results to save
ks =zeros(length(LAMV),length(NMV));
R2 = r2(NMV);
chis=zeros(length(LAMV),length(NMV));
d2gam2=zeros(length(LAMV),length(NMV));

%% start calculation
inm=1;
for NM=NMV
    ilam=1;
    for LAM=LAMV
        fprintf('NM=%d, LAM=%.2f // ', NM,LAM)
        [chis(ilam,inm),ks(ilam,inm),d2gam2(ilam,inm)]=spinodal(M,NM,LAM,FA);
        ilam=ilam+1;
    end
    inm=inm+1;
end

if NMV(1)>=1e4  % Gaussian chain limit
    save('data/spingc','LAMV','NMV','ks','chis','d2gam2')
elseif NMV(1)<=1e-3  % Rigid rod limit
else  % Worm-like chain
    save('data/spinwlc','LAMV','NMV','ks','chis','d2gam2')
end