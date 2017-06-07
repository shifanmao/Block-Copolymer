clear;
close all

% universal parameters
N = 100;

filename1 = '../data/figureS5A';
filename2 = '../data/figureS5B';
filename3 = '../data/figureS5C';
filename4 = '../data/figureS5D';
filename5 = '../data/figureS5E';
filename6 = '../data/figureS5F';
filenames = {filename1, filename2, filename3, filename4, filename5, filename6};

FAV = linspace(.1, .9, 11);
LAMV = linspace(-1, 1, 11);
CHIABSV = zeros(length(FAV), length(LAMV));
KSV = zeros(length(FAV), length(LAMV));

cnt = 1;
for PHIP = [.2, .5, .8]
    for NM = [1e-2, 1e2]
        for ii = 1:length(FAV)
            FA = FAV(ii)
            for jj = 1:length(LAMV)
                LAM = LAMV(jj);

                [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
                CHIABSV(ii, jj) = CHIABS*NM*PHIP;
                KSV(ii, jj) = KS*sqrt(r2(NM));
            end
        end

        filename = filenames{cnt};
        save(filename, 'CHIABSV', 'KSV', 'LAMV', 'FAV', 'PHIP', 'NM');
        cnt = cnt+1;
    end
end

% figure;surf(LAMV, FAV, KSV)


% 
% %%%% cuts %%%%
% clear;
% filename3 = '../data/figureS5C';
% load(filename3)
% 
% figure;hold
% for ii = 1:length(FAV)
%     plot(LAMV, KSV(ii,:))
% end
% 
% figure;hold
% for ii = 1:length(FAV)
%     plot(LAMV, CHIABSV(ii,:)*FAV(ii)*(1-FAV(ii))*4)
% end





% 
% 
% %%%% some conditions %%%%
% clear;
% N = 100;
% PHIP = .5;
% 
% FA = .2;
% NM = .01;
% % LAMV = linspace(-1, 0, 21);
% LAMV = linspace(-1, 1, 21);
% KSV = zeros(1, length(LAMV));
% CHIABSV = zeros(1, length(LAMV));
% 
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii)
%     [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
%     KSV(ii) = KS*sqrt(r2(NM));
%     CHIABSV(ii) = CHIABS*PHIP*NM;
% end
% 
% figure;plot(LAMV, CHIABSV)
% figure;plot(LAMV, KSV)
% 

%%%% structure factors %%%%
clear;
N = 100;
PHIP = .5;

FA = .2;
NM = .01;
LAM = -0.20;
KV = logspace(-3, 2, 100)/sqrt(r2(NM));

[CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
CHIV = linspace(.2, .9, 5);
KSV = zeros(length(CHIV), 1);

figure;hold
for ii = 1:length(CHIV)
    COL = (ii-1)/(length(CHIV)-1);
    
    CHI = CHIV(ii);
    CHIMAT = [0,CHIABS;CHIABS,0]*CHI
    
    [KS, EIGS, EIGVS] = eigs_wsolvent(N, NM, LAM, FA, PHIP, CHIMAT, 1);
    KSV(ii) = KS*sqrt(r2(NM));
    
    [EIG,EIGV]=gamma2_solvent(N,NM,LAM,FA,KV,CHIMAT,PHIP);
    plot(KV*sqrt(r2(NM)), 1/NM./EIG(:,1), '-', 'color', [COL 0 1-COL]);
%     plot(KV*sqrt(r2(NM)), 1/NM./EIG(:,2), '--', 'color', [COL 0 1-COL]);
end
set(gca,'xscale','log');set(gca,'yscale','log')