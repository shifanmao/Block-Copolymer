addpath(genpath('../chainstats/'))
clear;

% preset parametersI
N = 100;
NM = 100;
FA = 0.5;

NLAM = 51;
LAMV = linspace(-1,1,NLAM);
LAMV = 0;

% PHIPV = [0.5,1-1e-2];
PHIPV = linspace(0.1,1-1e-2,21);
NPHIP = length(PHIPV);

% results to return
CHIABS = zeros(NLAM,1);
CHIABKS = zeros(NLAM,1);
KSS = zeros(NLAM,1);
EIGS = zeros(NLAM,2);

% filename = 'solvent_spin.txt';
% fileID = fopen(filename,'wt');
% x = 0:.1:1;
% A = [x; exp(x)];
% 
% figure;hold
% for ii = 1:NPHIP
%     PHIP = PHIPV(ii)
%     
%     for jj = 1:length(LAMV)
%         LAM = LAMV(jj)
%         [CHIABS(jj),KSS(jj),EIGS(jj,:),CHIABKS(jj)] = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
%     end
%     pref = 4*PHIP*FA*(1-FA);
%     plot(LAMV(KSS*RM<=0.01),pref*CHIABS(KSS*RM<=0.01)*NM,'k--')
%     plot(LAMV(KSS*RM>0.01),pref*CHIABS(KSS*RM>0.01)*NM,'k-')
%     plot(LAMV,pref*CHIABKS*NM,'k-.')
% %     fprintf(,'%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f', LAM, FA, PHIP, CHIABS(jj),KSS(jj),EIGS(jj,:),CHIABKS(jj));
%     fprintf(fileID,'%.2f\n',CHIABS(jj));
% end
% fclose(fileID);




% write to file
filename='solvent_spin';
if ~exist(filename,'file')
    outfile = fopen(filename, 'wt');
    
    for NM = 100
        % WAVEVECTOR
        NK = 100;
        RM = (r2(NM))^0.5; % Normalization factor
        KV = logspace(-2,2,NK)/RM;  % Wavevector

        for FA = linspace(0.1,0.9,21)
            for ii = 1:NPHIP
                PHIP = PHIPV(ii);

                for jj = 1:length(LAMV)
                    LAM = LAMV(jj);
                    [CHIABS,KSS,EIGS,CHIABKS] = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
                    fprintf(outfile,'%.2f, %.4f, %.4f, %.4f, %.4f,     %.4f, %.4f, %.4f, %.4f\n', ...
                                     NM,   LAM,  FA,   PHIP, CHIABS*NM,KSS*RM,EIGS(1)*NM,EIGS(2)*NM,CHIABKS*NM);
                end
            end
        end
    end
    fclose(outfile);
end

% % make a phase diagram of RBC in solvent
% figure;hold
% plot(PHIPV(KSS*RM<=0.01),CHIABS(KSS*RM<=0.01)*NM,'k--')
% plot(PHIPV(KSS*RM>0.01),CHIABS(KSS*RM>0.01)*NM,'k-')














NK = 100;
RM = (r2(NM))^0.5; % Normalization factor
KV = logspace(-2,2,NK)/RM;  % Wavevector
LAM = 0;


%%% stop

NPHIP = 18;
PHIPV = linspace(0.01,0.99,NPHIP);

NFA = 16;
FAV = linspace(0.01,0.99,NFA);

% CHIABS = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
CHIAB = 0/NM; % Flory-Huggins factor between A and S

X1 = zeros(length(FAV), length(PHIPV));
X2 = zeros(length(FAV), length(PHIPV));
V1 = zeros(length(FAV), length(PHIPV));
V2 = zeros(length(FAV), length(PHIPV));
PREF = zeros(length(FAV), length(PHIPV));
CHIABS = zeros(length(FAV), length(PHIPV));
for ii = 1:length(FAV)
    FA = FAV(ii)
    for jj = 1:length(PHIPV)
        PHIP = PHIPV(jj);
        
        % CHI PARAMETERS
        CHIBS = 0; % Flory-Huggins factor between B and S
        CHIAS = 0;
        CHIBA = CHIAB;
        CHI = [CHIAS,CHIAB;CHIBA,CHIBS];

        [EIG1,EIG2,EIGV1,EIGV2,KS1,KS2]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,PHIP);
        KSIND = find(KV >= KS1, 1);
        if EIGV1(KSIND,1)+EIGV1(KSIND,2)<0
            EIGV1(KSIND,:) = -EIGV1(KSIND,:);
        end
        
        X1(ii,jj) = PHIP;
        X2(ii,jj) = FA*PHIP;
        V1(ii,jj) = (EIGV1(KSIND,1)+EIGV1(KSIND,2));
        V2(ii,jj) = EIGV1(KSIND,1);
        PREF(ii,jj) = (1/EIG1(KSIND)/NM);
        if PREF(ii,jj) > 100
            PREF(ii,jj) = nan;
        end
        CHIABS(ii,jj) = solvnt_spin(N,NM,LAM,FA,KV,PHIP);
    end
end

figure;hold
hq = quiver(X1,X2,V1.*PREF,V2.*PREF,'color','k', 'linewidth', 1.);
ABOVESPIN = find(PREF == nan);
plot(X1(ABOVESPIN), X2(ABOVESPIN), 'r.', 'markersize', 20)

% headWidth = 8;
% headLength = 8;
% LineLength = 0.08;
% 
% %get the data from regular quiver
% U = hq.UData;
% V = hq.VData;
% X = hq.XData;
% Y = hq.YData;
% 
% for ii = 1:length(X)
%     for ij = 1:length(X)
% 
%         headWidth = 5;
%         ah = annotation('arrow',...
%             'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
%         set(ah,'parent',gca);
%         set(ah,'position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);
% 
%     end
% end

plot([0,1], [0,1], 'k-', 'linewidth', 1.5)
plot([0,1], [0,0.5], 'k-' , 'linewidth', 1.5)
axis([0,1,0,1])
set(gca,'fontsize',30)
xlabel('\phi_P')
ylabel('\phi_A=\phi_Pf_A')
box on

figure;surf(X1,X2,CHIABS*NM)


figure;hold
for I = 2:length(PHIPV)
    COL = (I-1)/(length(PHIPV)-1);
    plot(FAV, CHIABS(:,I)*NM,'-', 'linewidth', 2, ...
        'color', [COL 0 1-COL])
end