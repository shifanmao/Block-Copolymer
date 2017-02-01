function [EIG1,EIG2,EIGV1,EIGV2,THETA1,THETA2,KS1,KS2]=gamma2_solvent(N,NM,LAM,FA,kv,CHI,FP)

%% Inputs
FS = 1-FP;   % Fraction of solvent

%% Start Calculation
%%% CALCULATE S FUNCTIONS %%%
EIG1 = zeros(length(kv),1);
EIG2 = zeros(length(kv),1);
EIGV1 = zeros(length(kv),2);
EIGV2 = zeros(length(kv),2);

for ii = 1:length(kv)
  k = kv(ii);  % wavevector
  s2inv = s2inverse(N,NM,LAM,FA,k);
  SAAINV = s2inv(1,1);
  SABINV = s2inv(1,2);
  SBBINV = s2inv(2,2);
  
  %%% CALCULATE GAM FUNCTIONS %%%
  GAMAA = -2*CHI(1,1)+N*NM*SAAINV./FP+1/FS;
  GAMBB = -2*CHI(2,2)+N*NM*SBBINV./FP+1/FS;
  GAMAB = (CHI(1,2)-CHI(1,1)-CHI(2,2))+N*NM*SABINV/FP+1/FS;
  GAMBA = GAMAB;
  GAM = [GAMAA, GAMAB; GAMBA, GAMBB];

  %%% CALCULATE EIGENMODES %%%
  [EIV,EI] = eig(GAM);
  EI = [EI(1,1),EI(2,2)];
  [~,ind] = sort(EI,'ascend');ind1 = ind(1);ind2 = ind(2);
  EIG1(ii) = EI(ind1);
  EIG2(ii) = EI(ind2);
  EIGV1(ii,:) = EIV(:, ind1);
  EIGV2(ii,:) = EIV(:, ind2);
end

THETA1 = atan(EIGV1(1,1)/EIGV1(1,2));
THETA2 = atan(EIGV2(1,1)/EIGV2(1,2));

%%% FIND CRITICAL WAVEMODE %%%
KS1 = kv(EIG1==min(EIG1));
KS2 = kv(EIG2==min(EIG2));
end
