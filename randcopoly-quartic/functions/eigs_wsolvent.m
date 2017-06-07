function [KS, EIGS, EIGVS] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHI, IEIG)

RM = sqrt(r2(NM));

function EIGMIN = eig_wsolvent(K)
   [vals,~]=gamma2_solvent(N,NM,LAM,FA,K,CHI,PHIP);
   EIGMIN = vals(:,IEIG);
end

Kmin = 1e-2/RM;
Kmax = 1e1/RM;
% % take the best of two minimizers
KS1 = fminbnd(@eig_wsolvent,Kmin,Kmax);  % can get stuck at finite k
KS2 = Kmin;
% KS2 = fminsearch(@eig_wsolvent,1/RM);    % can get stuck at k=0
if eig_wsolvent(KS1) < eig_wsolvent(KS2)
    KS = KS1;
else
    KS = KS2;
end

[EIG,EIGV]=gamma2_solvent(N,NM,LAM,FA,KS,CHI,PHIP);
EIGVS = EIGV(IEIG*2 - 1:IEIG*2);
EIGS = EIG(IEIG);

end