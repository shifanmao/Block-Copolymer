function CHIABS = solvnt_spin(N,NM,LAM,FA,KV,FP)

chis = 0/NM;
chie = 100/NM;
while (chie - chis) > 1e-2/NM
    mid = (chie + chis) / 2;
    CHIAB = mid;
    CHIBA = CHIAB;
    CHI = [0, CHIAB; CHIBA, 0];

    [EIG1, ~, ~, ~, KS1, ~] = gamma2_solvent(N,NM,LAM,FA,KV,CHI,FP);
    KSIND = find(KV >= KS1, 1);
    if EIG1(KSIND)*NM < 1e-2
        chie = mid;
    else
        chis = mid;
    end
end
CHIABS = (chis + chie) / 2;

end
