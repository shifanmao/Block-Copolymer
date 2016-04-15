% This plots critical wavemode and peak widths
% from mean-field theories and Monte-Carlo simulations
clear;close all

% start code
LAMV_SIM = -0.75:0.25:0.25;
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
Lbox=20;
Lbin=1;
ICHI=41;

figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',20)

% plot mean-field theory results
ieps = 1;
for EPS = EPSV
    col = (ieps-1)/(length(EPSV)-1);
    
    NM=G*EPS;  % number of Kuhn steps per monomer
    data = load(sprintf('data/WLC_NM%.2f',NM));
    LAMV_MF = data(:,1);
    KS_MF = data(:,2);
    ALPHA_MF = data(:,4);
    
    LIFT=find(KS_MF>1e-3);LIFT=LIFT(end);
    figure(1);plot(LAMV_MF,KS_MF,'--','linewidth',3,'color',[col 0 1-col])
%               plot(LAMV_MF(LIFT),0,'o','MarkerEdgeColor',[col 0 1-col],'MarkerFaceColor',[col 0 1-col],...
%                   'markersize',12,'color',[col 0 1-col],'linewidth',2);
    figure(2);plot(LAMV_MF,ALPHA_MF,'--','linewidth',3,'color',[col 0 1-col])
    ieps = ieps+1;
end

% plot simulation results
ieps = 1;p1 = zeros(length(EPSV),1);p2 = zeros(length(EPSV),1);
for EPS = EPSV
    col = (ieps-1)/(length(EPSV)-1);
    
    KSV_SIM = zeros(length(LAMV_SIM),46);
    ALPHA_SIM = zeros(length(LAMV_SIM),46);
    ERRV_SIM = zeros(length(LAMV_SIM),46,3);
    ilam=0;
    for LAM = LAMV_SIM
        ilam = ilam+1;
        [~,~,~,KS_SIM,~,~,~,D2S_SIM,ERR_SIM]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        KSV_SIM(ilam,:) = KS_SIM;
        ALPHA_SIM(ilam,:) = sqrt(D2S_SIM/2);
        ERRV_SIM(ilam,:,1) = ERR_SIM(:,1);
        ERRV_SIM(ilam,:,2) = ERR_SIM(:,2);
        ERRV_SIM(ilam,:,3) = ERR_SIM(:,3);
    end
    
%    figure(1);p1(ieps)=plot(LAMV_SIM,KSV_SIM(:,ICHI),'o',...
%        'MarkerEdgeColor',[col 0 1-col],'markersize',12,'linewidth',2);
%    figure(2);p2(ieps)=plot(LAMV_SIM,ALPHA_SIM(:,ICHI),'o',...
%        'MarkerEdgeColor',[col 0 1-col],'markersize',12,'linewidth',2);

    figure(1);p1(ieps)=errorbar(LAMV_SIM,KSV_SIM(:,ICHI),ERRV_SIM(:,ICHI,1),'o',...
	'MarkerEdgeColor',[col 0 1-col],'markersize',12,'color',[col 0 1-col],'linewidth',2);
    figure(2);p2(ieps)=errorbar(LAMV_SIM,ALPHA_SIM(:,ICHI),ERRV_SIM(:,ICHI,3),'o',...
	'MarkerEdgeColor',[col 0 1-col],'markersize',12,'color',[col 0 1-col],'linewidth',2);
    ieps = ieps+1;
end

figure(1);
xlim([-1,.5]);ylim([0,5])
xlabel('\lambda');ylabel('Critical Wavemode R_Mq^*');box on
set(gca,'Xtick',-1:0.25:0.5)
set(gca,'XtickLabel',{'-1','-0.75','-0.5','-0.25','0','0.25','0.5'})
l1=legend(p1,{'N_M=0.05','N_M=0.50','N_M=5.00'},'location','northeast');
M = findobj(l1,'type','line');
set(M,'linewidth',2)

figure(2);
xlabel('\lambda');
ylabel('Inverse Peak Width \alpha');
box on
xlim([-1,.5]);
ylim([0,3])
set(gca,'yscale','linear')
% set(gca,'yscale','log');
% ylim([1e-2,1e2])
set(gca,'Xtick',-1:0.25:0.5)
set(gca,'XtickLabel',{'-1','-0.75','-0.5','-0.25','0','0.25','0.5'})
l2=legend(p2,{'N_M=0.05','N_M=0.50','N_M=5.00'},'location','northeast');
M = findobj(l2,'type','line');
set(M,'linewidth',2)

% end code

figure(1);
savename = sprintf('../../results/randcopoly-results/random-simulation/figure5A.eps');
saveas(gcf,savename,'epsc')

figure(2);
savename = sprintf('../../results/randcopoly-results/random-simulation/figure5B.eps');
saveas(gcf,savename,'epsc')
