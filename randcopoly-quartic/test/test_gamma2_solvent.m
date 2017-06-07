clear;

% universal parameters
N = 100;
NM = 1e-2;

% LAM = -0.75;
% CHI = 0;

% [FAVV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, NM, LAM, CHI, 1, 0);
% figure;plotphase_wsolvent(FAVV, PHIPV, EIGV, EIG, KSV, NM)


CHI = [0,0;0,0];
LAM = -0.75;
FA = 0.2;
KV = logspace(-2, 5, 100);
PHIP = .5;

[CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
CHIABS = 0;
CHIMAT = [0,CHIABS;CHIABS,0];

[EIG,EIGV]=gamma2_solvent(N,NM,LAM,FA,KV,CHIMAT,PHIP);
figure;loglog(KV, 1./EIG(:,1)/NM);

% FA = 0.5;
% PHIP = 0.5;
% 
% LAMV = linspace(-1, 1, 20);
% EIGSV = zeros(length(LAMV), 1);
% CHIABSV = zeros(length(LAMV), 1);
% KSV = zeros(length(LAMV), 1);
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii)
% %     [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
% %     EIGSV(ii) = EIGS;
% %     CHIABSV(ii) = CHIABS;
%     
%     CHI = [0, 0;0, 0];
%     [KS, EIGS, EIGVS] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHI, 1);
%     EIGSV(ii) = EIGS;
%     KSV(ii) = KS;
% end


% % critical wavemode
% LAMV = linspace(-1, 1, 100);
% KSV = zeros(1, length(LAMV));
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii);
%     [KSV(ii), EIGS, EIGVS] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHI, 1);
% end
% 
% % spinodals
% FAV = linspace(0.01, .99, 100);
% CHIABS = zeros(1, length(FAV));
% KS = zeros(1, length(FAV));
% PHIP = 0.2;
% 
% LAM = -1;
% for ii = 1:length(FAV)
%     FA = FAV(ii);
%     [CHIABS(ii), KS(ii), EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
% end
% figure;plot(FAV, CHIABS*NM*PHIP);ylim([0,50])
% figure;plot(FAV, KS*sqrt(r2(NM)))

% 
% % spinodals
% LAMV = linspace(-1, 1, 51);
% CHIABS = zeros(1, length(LAMV));
% KS = zeros(1, length(LAMV));
% FA = 0.5;
% PHIP = 1;
% 
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii)
%     [CHIABS(ii), KS(ii), EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
% end
% figure;plot(LAMV, CHIABS*NM*PHIP);
% figure;plot(LAMV, KS*sqrt(r2(NM)))
% 



% % spinodals
% N = 100;
% 
% FA = .5;
% for NM = [1, 10, 100, 1e3]
%     figure;hold
%     for PHIP = [.5, 1.0]
%         disp(sprintf('NM=%d, FA=%.2f, PHIP=%.2f', NM, FA, PHIP))
%         LAMV = linspace(-1, 1, 51);
%         CHIABS = zeros(1, length(LAMV));
%         KS = zeros(1, length(LAMV));
% 
%         for ii = 1:length(LAMV)
%             LAM = LAMV(ii);
%             [CHIABS(ii), KS(ii), EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
%         end
%         
%         if PHIP == .5
%             plot(LAMV, CHIABS*NM*PHIP*4*FA*(1-FA), LAMV, KS*sqrt(r2(NM)));
%         else
%             plot(LAMV, CHIABS*NM*PHIP*4*FA*(1-FA), '--', LAMV, KS*sqrt(r2(NM)), '--');
%         end
%         title(sprintf('NM=%d, FA=%.2f, PHIP=%.2f', NM, FA, PHIP))
%     end
% end
% 
% 
% 
% 
% 
% 
% % spinodals
% N = 100;
% 
% for PHIP = [.5, 1.0];
% for NM = [1, 10, 100, 1e3]
%     figure;hold
%     for FA = [.2, .5]
%         disp(sprintf('NM=%d, FA=%.2f, PHIP=%.2f', NM, FA, PHIP))
%         LAMV = linspace(-1, 1, 51);
%         CHIABS = zeros(1, length(LAMV));
%         KS = zeros(1, length(LAMV));
% 
%         for ii = 1:length(LAMV)
%             LAM = LAMV(ii);
%             [CHIABS(ii), KS(ii), EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
%         end
%         
%         if FA == .5
%             plot(LAMV, CHIABS*NM*PHIP*4*FA*(1-FA), LAMV, KS*sqrt(r2(NM)));
%         else
%             plot(LAMV, CHIABS*NM*PHIP*4*FA*(1-FA), '--', LAMV, KS*sqrt(r2(NM)), '--');
%         end
%         title(sprintf('NM=%d, FA=%.2f, PHIP=%.2f', NM, FA, PHIP))
%     end
% end
% end


% % spinodals
% FAV = linspace(0.01, .99, 51);
% LAMV = linspace(-1, 1, 51);
% CHIABS = zeros(length(LAMV), length(FAV));
% KS = zeros(length(LAMV), length(FAV));
% PHIP = 0.2;
% 
% LAM = -1;
% for ii = 1:length(FAV)
%     FA = FAV(ii);
%     [CHIABS(ii), KS(ii), EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
% end
% figure;plot(FAV, CHIABS*NM*PHIP);ylim([0,50])
% figure;plot(FAV, KS*sqrt(r2(NM)))