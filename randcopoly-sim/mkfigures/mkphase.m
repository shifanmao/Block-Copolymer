% make phase diagrams with MF predictions and MC simulations
clear;
close all;

addpath('../functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%% SIM # 1 : NM = 0.05 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% NM=0.05; % each monomer has 0.1 Kuhn steps
for NM = [0.05,0.50,5.00]
COL='b';

N=8;  % total of 8 monomers
G=5;
FA=0.5;    % equal chemical composition
RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers

% find spinodal CHIS
LAMV=linspace(-1,.99,201); % anti-correlated random copolymer
KS = zeros(length(LAMV),1);
CHIS = zeros(length(LAMV),1);

for ii = 1:length(LAMV)
    LAM = LAMV(ii);
    fprintf('NM=%.2f, LAM=%.2f\n',NM, LAM) 
    [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
    KS(ii)=kval;
    CHIS(ii)=0.5*sval;  % spinodal
    D2S(ii)=1/(sval^2*RM^2)*d2gam2;
end
% 
% figure;plot(LAMV,RM*KS,'-','color',COL)
% xlabel('\lambda');ylabel('R_Mq^*');box on

% find Lifshitz point
KSP = RM*KS;
IND = find(KSP<=0.01,1);
LAML = LAMV(IND)

result = [LAML,0];
result = [result; [LAMV',CHIS*NM]];

filename = strcat('simdata/',...
    sprintf('MFspinNM%.2f',NM));
dlmwrite(filename,result);
end



% figure;hold;set(gca,'fontsize',20)
% plot(LAMV,CHIS*NM,'--','color',COL,'linewidth',2)
% plot([LAML,LAML],[0,10],'-.','color',COL,'linewidth',2)
% xlabel('\lambda');ylabel('\chi_SvN_M');
% xlim([-1,.50]);box on
% 
% % load simulation results
% x = load('simdata/phasediag');
% INDSIM = find(x(:,1)==NM/G);
% LAMODT  = x(INDSIM,2)';
% CHI0 = zeros(length(LAMODT),1)';
% CHIM = ones(length(LAMODT),1)'*10;
% CHIODT1 = x(INDSIM,3)';
% CHIODT2 = x(INDSIM,4)';
% 
% plot(LAMODT,CHIODT1,'b.-','markersize',20)
% plot(LAMODT,CHIODT2,'r.-','markersize',20)
% 
% % fill areas
% h1 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHI0(3:end),fliplr(CHIODT1(3:end))],'k');
% h2 = fill([LAMODT(1:3),fliplr(LAMODT(1:3))],[CHI0(1:3),fliplr(CHIODT2(1:3))],'k');
% alpha(.2);
% set(h1,'EdgeColor','None');
% set(h2,'EdgeColor','None');
% 
% h3 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHIODT2(3:end),fliplr(CHIODT1(3:end))],...
%      'b');
% alpha(.2);set(h3,'EdgeColor','None');
% 
% h3 = fill([LAMODT,fliplr(LAMODT)],[CHIODT2,fliplr(CHIM)],...
%      'r');
% alpha(.2);set(h3,'EdgeColor','None');



% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%% SIM # 1 : NM = 0.50 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% NM=0.50; % each monomer has 0.1 Kuhn steps
% COL='b';
% 
% N=8;  % total of 8 monomers
% G=5;
% FA=0.5;    % equal chemical composition
% RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
% 
% % find spinodal CHIS
% LAMV=linspace(-1,.99,201); % anti-correlated random copolymer
% KS = zeros(length(LAMV),1);
% CHIS = zeros(length(LAMV),1);
% 
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii);
%     fprintf('NM=%.2f, LAM=%.2f\n',NM, LAM) 
%     [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
%     KS(ii)=kval;
%     CHIS(ii)=0.5*sval;  % spinodal
%     D2S(ii)=1/(sval^2*RM^2)*d2gam2;
% end
% 
% figure;plot(LAMV,RM*KS,'-','color',COL)
% xlabel('\lambda');ylabel('R_Mq^*');box on
% 
% % find Lifshitz point
% KSP = RM*KS;
% IND = find(KSP<=0.01,1);
% LAML = LAMV(IND)
% 
% figure;hold;set(gca,'fontsize',20)
% plot(LAMV,CHIS*NM,'--','color',COL,'linewidth',2)
% plot([LAML,LAML],[0,10],'-.','color',COL,'linewidth',2)
% xlabel('\lambda');ylabel('\chi_SvN_M');
% xlim([-1,.50]);box on
% 
% % load simulation results
% x = load('simdata/phasediag');
% INDSIM = find(x(:,1)==NM/G);
% LAMODT  = x(INDSIM,2)';
% CHI0 = zeros(length(LAMODT),1)';
% CHIM = ones(length(LAMODT),1)'*10;
% CHIODT1 = x(INDSIM,3)';
% CHIODT2 = x(INDSIM,4)';
% 
% plot(LAMODT,CHIODT1,'b.-','markersize',20)
% plot(LAMODT,CHIODT2,'r.-','markersize',20)
% 
% % fill areas
% h1 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHI0(3:end),fliplr(CHIODT1(3:end))],'k');
% h2 = fill([LAMODT(1:3),fliplr(LAMODT(1:3))],[CHI0(1:3),fliplr(CHIODT2(1:3))],'k');
% alpha(.2);
% set(h1,'EdgeColor','None');
% set(h2,'EdgeColor','None');
% 
% h3 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHIODT2(3:end),fliplr(CHIODT1(3:end))],...
%      'b');
% alpha(.2);set(h3,'EdgeColor','None');
% 
% h3 = fill([LAMODT,fliplr(LAMODT)],[CHIODT2,fliplr(CHIM)],...
%      'r');
% alpha(.2);set(h3,'EdgeColor','None');
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%% SIM # 1 : NM = 5.00 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% NM=5.00; % each monomer has 0.1 Kuhn steps
% COL='b';
% 
% N=8;  % total of 8 monomers
% G=5;
% FA=0.5;    % equal chemical composition
% RM=sqrt(r2wlc(NM));  % end-to-end distance of a monomers
% 
% % find spinodal CHIS
% LAMV=linspace(-1,.99,201); % anti-correlated random copolymer
% KS = zeros(length(LAMV),1);
% CHIS = zeros(length(LAMV),1);
% 
% for ii = 1:length(LAMV)
%     LAM = LAMV(ii);
%     fprintf('NM=%.2f, LAM=%.2f\n',NM, LAM) 
%     [kval,sval,d2gam2]=kmaxwlc(N,NM,FA,LAM);
%     KS(ii)=kval;
%     CHIS(ii)=0.5*sval;  % spinodal
%     D2S(ii)=1/(sval^2*RM^2)*d2gam2;
% end
% 
% figure;plot(LAMV,RM*KS,'-','color',COL)
% xlabel('\lambda');ylabel('R_Mq^*');box on
% 
% % find Lifshitz point
% KSP = RM*KS;
% IND = find(KSP<=0.01,1);
% LAML = LAMV(IND)
% 
% figure;hold;set(gca,'fontsize',20)
% plot(LAMV,CHIS*NM,'--','color',COL,'linewidth',2)
% plot([LAML,LAML],[0,10],'-.','color',COL,'linewidth',2)
% xlabel('\lambda');ylabel('\chi_SvN_M');
% xlim([-1,.50]);box on
% 
% % load simulation results
% x = load('simdata/phasediag');
% INDSIM = find(x(:,1)==NM/G);
% LAMODT  = x(INDSIM,2)';
% CHI0 = zeros(length(LAMODT),1)';
% CHIM = ones(length(LAMODT),1)'*10;
% CHIODT1 = x(INDSIM,3)';
% CHIODT2 = x(INDSIM,4)';
% 
% plot(LAMODT,CHIODT1,'b.-','markersize',20)
% plot(LAMODT,CHIODT2,'r.-','markersize',20)
% 
% % fill areas
% h1 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHI0(3:end),fliplr(CHIODT1(3:end))],'k');
% h2 = fill([LAMODT(1:3),fliplr(LAMODT(1:3))],[CHI0(1:3),fliplr(CHIODT2(1:3))],'k');
% alpha(.2);
% set(h1,'EdgeColor','None');
% set(h2,'EdgeColor','None');
% 
% h3 = fill([LAMODT(3:end),fliplr(LAMODT(3:end))],[CHIODT2(3:end),fliplr(CHIODT1(3:end))],...
%      'b');
% alpha(.2);set(h3,'EdgeColor','None');
% 
% h3 = fill([LAMODT,fliplr(LAMODT)],[CHIODT2,fliplr(CHIM)],...
%      'r');
% alpha(.2);set(h3,'EdgeColor','None');
