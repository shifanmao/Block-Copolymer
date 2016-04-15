addpath('../functions')

N=100;  % total of 100 monomers
FA=0.5;    % equal chemical composition

% Gaussian chain
NM=100; % each monomer has 100 Kuhn steps
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

LAMV = linspace(-.99,.99,501);
KS = zeros(length(LAMV),1);
CHIS = zeros(length(LAMV),1);
ALPHA = zeros(length(LAMV),1);
for ii = 1:length(LAMV)
    LAM = LAMV(ii);
    [kval,sval,d2gam2]=kmaxgc(N,NM,FA,LAM);
    KS(ii)=kval;
    CHIS(ii)=0.5*sval;  % spinodal
    ALPHA(ii)=d2gam2./(RM^2);
end

data = [LAMV',RM*KS,CHIS*NM,sqrt(ALPHA*NM/2)];
dlmwrite(sprintf('data/GC'),data)

% WLC
NMV = logspace(-2,2,11);
LAML_WLC = zeros(length(NMV),1);
cnt = 1;
NMV = [0.05,0.5];

for NM = NMV
    col = (cnt-1)/(length(NMV)-1);
    NM
    RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

    LAMV = linspace(-.99,.99,501);
    KS = zeros(length(LAMV),1);
    CHIS = zeros(length(LAMV),1);
    ALPHA = zeros(length(LAMV),1);
    for ii = 1:length(LAMV)
        LAM = LAMV(ii)
        [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
        KS(ii)=kval;
        CHIS(ii)=0.5*sval;  % spinodal
        ALPHA(ii)=d2gam2./(RM^2);
    end
    
    data = [LAMV',RM*KS,CHIS*NM,sqrt(ALPHA*NM/2)];
    dlmwrite(sprintf('data/WLC_NM%.2f',NM),data)
    cnt = cnt+1;
end
% 
% % Rigid rod
% N = 100;
% NM=0.01; % each monomer has 100 Kuhn steps
% RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
% 
% LAMV = linspace(-.99,.99,501);
% KS = zeros(length(LAMV),1);
% CHIS = zeros(length(LAMV),1);
% ALPHA = zeros(length(LAMV),1);
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii)
%     [kval,sval,d2gam2]=kmaxrr(N,NM,FA,LAM);
%     KS(ii)=kval;
%     CHIS(ii)=0.5*sval;  % spinodal
%     ALPHA(ii)=d2gam2./(RM^2);
% end
% 
% data = [LAMV',RM*KS,CHIS*NM,sqrt(ALPHA*NM/2)];
% dlmwrite(sprintf('data/RR'),data)