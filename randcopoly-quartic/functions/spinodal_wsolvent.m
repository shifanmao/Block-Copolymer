function [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP)

chis = 0/NM/PHIP;
chie = 3e1/NM/PHIP;

while (chie - chis) > 1e-2/NM/PHIP
    mid = (chie + chis) / 2;
    CHI = setCHIAB(mid);

    [~, MINEIG, ~] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHI, 1);
    if MINEIG*N < 0
        chie = mid;
    else
        chis = mid;
    end
end
CHIABS = chis;

CHI = setCHIAB(CHIABS);
[KS, EIGS, EIGVS] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHI, 1);

% reset KS
if KS*sqrt(r2(NM)) <= 1e-2
    KS = 0;
end

end

function CHI = setCHIAB(CHIAB)
    CHIBA = CHIAB;
    CHI = [0, CHIAB; CHIBA, 0];
end
