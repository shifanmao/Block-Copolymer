function R2=r2(NM)
% calculates average end-to-end squared distance
% of polymers with NM number of Kuhn steps

if NM>=1e4  % Gaussian chain limit
    R2 = NM;
elseif NM<=1e-4  % Rigid rod limit
    R2 = NM^2;
else  % Worm-like chain
    R2 = NM-(1/2)*(1-exp(-2*NM));
end