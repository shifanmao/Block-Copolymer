% This plots radius of gyrations calculated from
% mean-field theory and Monte-Carlo simulation
clear;close all

for LAM=[-0.75,0.00]
    EPSV = [0.01,0.10,1.00];
    PLOTON = 1;

    % simulation constants
    M=8;        % number of blocks
    G=5;  % number of discrete monomers
    N=M*G;

    %f=figure('Position', [100, 100, 900, 900]);
    figure;
    hold;set(gca,'fontsize',18)
    NK=linspace(0,8,1000);

    for EPS = EPSV
        plotsim3(EPS,LAM,PLOTON);
    end
    
    for EPS = EPSV
        EPSV = [0.01,0.10,1.00];
        coli = find(EPSV==EPS);
        col = (coli-1)/(length(EPSV)-1);
        
        NM = EPS*G;
        DM = EPS*NK*G;
        RM = NM-(1/2)*(1-exp(-2*NM));
        RD = DM-(1/2)*(1-exp(-2*DM));
        
%         plot(NK,sqrt(RD/RM),'k-','linewidth',2)
        plot(NK,sqrt(RD/RM),'-','color',[col 0 1-col],'linewidth',2)
    end
    
    plot([1,1],[0.1,1],'k-','linewidth',2)
    plot([0.1,1],[1,1],'k-','linewidth',2)
    
    xlabel('J/G');ylabel('R_{J}/R_{M}')
    box on
    
    set(gca,'Xtick',0:8)
    set(gca,'XtickLabel',{'0','1','2','3','4','5','6','7','8'})

    % add power laws
    xlog = linspace(5,8,5);
    ylog = power(xlog,1)*1.1;
    plot(xlog,ylog,'k--','linewidth',2)
    
    xlog = linspace(5,8.6,5);
    ylog = power(xlog,1/2)*0.95;
    plot(xlog,ylog,'k--','linewidth',2)
    text(5.5,8,'N','FontSize',20)
    text(6,2,'N^{1/2}','FontSize',20)
    
    xlim([0.1,10]);ylim([0.1,10])

    if LAM==-0.75
        savename = sprintf('../../results/randcopoly-results/random-simulation/figure7A.eps');
    elseif LAM==0.00
        savename = sprintf('../../results/randcopoly-results/random-simulation/figure7B.eps');
    else
        error('LAM undefined')
    end
    
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    legend_str=[{'N_M=0.05 (\chi=4\chi_{S}^{MF})'},{'N_M=0.05 (\chi=0)'},...
                {'N_M=0.50 (\chi=4\chi_{S}^{MF})'},{'N_M=0.50 (\chi=0)'},...
                {'N_M=5.00 (\chi=4\chi_{S}^{MF})'},{'N_M=5.00 (\chi=0)'}];
    %legend(legend_str,'location','northwest')
    saveas(gcf,savename,'epsc')
end