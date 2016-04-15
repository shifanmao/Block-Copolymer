% Plot SINV from mean-field theory and Monte-Carlo simulations
clear;close all

% start code
LAMV = [-0.75,0.00];
EPSV = [0.01,0.10,1.00];
PLOTON = 0;

% simulation parameters
N=8;
G=5;
RM=2;
Lbox=20;
Lbin=1;

% plot range
ind = 1:2:41;

for iLAM=1:length(LAMV)
    LAM=LAMV(iLAM);
    figure(iLAM);hold;set(gca,'fontsize',20)

    % plot mean-field theory results
    ieps = 1;
    for EPS = EPSV
        col = (ieps-1)/(length(EPSV)-1);
        
        [chis,chiv,ks,ksim,SINV_MF,SINV_SIM,ERR]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        
        % plot range
        chiv = chiv(ind);
        SINV_SIM = SINV_SIM(ind);
        SINV_MF = SINV_MF(ind);
        if LAM>=0
            plot(chiv*G,SINV_MF,'--','linewidth',3,'color',[0 0 0])
        else
            plot(chiv*G,SINV_MF,'--','linewidth',3,'color',[col 0 1-col])
        end
        ieps = ieps+1;
    end

    % plot simulation results
    cnt=1;p1=zeros(length(EPSV),1);
    ieps = 1;
    for EPS = EPSV
        col = (ieps-1)/(length(EPSV)-1);
        
        [chis,chiv,ks,ksim,SINV_MF,SINV_SIM,~,~,ERR]=plotsim(EPS,LAM,PLOTON);
        NM=EPS*G;
        R2=-0.5+0.5*exp(-2*NM)+NM;
        
        % plot range
        chiv = chiv(ind);
        SINV_SIM = SINV_SIM(ind);
        ksim = ksim(ind);
        err = ERR(ind,2);
%        p1(cnt)=plot(chiv*G,SINV_SIM,'o',...
%            'MarkerEdgeColor',[col 0 1-col],'linewidth',3,'markersize',10,'linewidth',2);
        p1(cnt)=errorbar(chiv*G,SINV_SIM,err,'o',...
            'MarkerEdgeColor',[col 0 1-col],'color',[col 0 1-col],...
	    'linewidth',3,'markersize',10,'linewidth',2);
        ieps = ieps+1;
        cnt = cnt+1;
    end

    xlabel('\chivG');ylabel('S^{-1}(q^*)');box on
    l1=legend(p1,{'N_M=0.05','N_M=0.50','N_M=5.00'},'location','northeast');
        %'Position',[0.72,0.25,0.1,0.1])
    M = findobj(l1,'type','line');
    set(M,'linewidth',2)
    if LAM==0
        ylim([0,1.2]);xlim([-0.1,8]);
        set(gca,'Ytick',0:0.2:1.2)
        set(gca,'YtickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2'})
    elseif LAM==-0.75
        ylim([0,2.4]);xlim([-0.1,20]);
        set(gca,'Ytick',0:0.4:2.4)
        set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0','2.4'})
    else
        ylim([0,2]);xlim([0,20]);
        set(gca,'Ytick',0:0.4:2.0)
        set(gca,'YtickLabel',{'0.0','0.4','0.8','1.2','1.6','2.0'})
    end

    % end code
    %savename =
    %sprintf('../../results/randcopoly-results/random-simulation/ssim-lam%.2f.eps',LAM);
end

figure(1);
savename = '../../results/randcopoly-results/random-simulation/figure6A.eps';
saveas(gcf,savename,'epsc')

figure(2);
savename = '../../results/randcopoly-results/random-simulation/figure6B.eps';
saveas(gcf,savename,'epsc')