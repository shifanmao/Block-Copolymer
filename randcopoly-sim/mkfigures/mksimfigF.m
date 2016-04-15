% This plots structure factors calculated from
% mean-field theory and Monte-Carlo simulation
clear;close all

LAMV = -0.75;
EPSV = [0.01,0.10,1.00];
PLOTON = 1;

figure(1);hold;set(gca,'fontsize',18)
figure(2);hold;set(gca,'fontsize',18)
NK=logspace(-2,2,20);
R2 = NK-(1/2)*(1-exp(-2*NK));
plot(NK,sqrt(R2),'k-','linewidth',2)

for LAM = LAMV
    for EPS = EPSV
        plotsim3(EPS,LAM,PLOTON);
    end
end

figure(2);
ylim([1e-2,1e1])

% add power laws
x = logspace(-2,-1.5,5);y = 1.5*x.^1;
plot(x,y,'k--','linewidth',2)
x = logspace(1.3,2,5);y = 0.7*x.^(1/2);
plot(x,y,'k--','linewidth',2)

xlabel('\epsilonN');
ylabel('R_{N}')
% ylabel('$<|\vec{r}_{iN}-\vec{r}_{i1}|^2>$','Interpreter','latex','FontName','Arial')
set(gca,'xscale','log');set(gca,'yscale','log');
box on
text(1.2e-2,3.5e-2,'N','FontSize',20)
text(40,3.5,'N^{1/2}','FontSize',20)

savename = sprintf('../../results/randcopoly-results/random-simulation/figure7A.eps');
saveas(gcf,savename,'epsc')


% This plots structure factors calculated from
% mean-field theory and Monte-Carlo simulation
clear;close all

LAMV = 0;
EPSV = [0.01,0.10,1.00];
PLOTON = 1;

figure(1);hold;set(gca,'fontsize',18)
figure(2);hold;set(gca,'fontsize',18)
NK=logspace(-2,2,20);
R2 = NK-(1/2)*(1-exp(-2*NK));
plot(NK,sqrt(R2),'k-','linewidth',2)

for LAM = LAMV
    for EPS = EPSV
        plotsim3(EPS,LAM,PLOTON);
    end
end

figure(2);
ylim([1e-2,1e1])

% add power laws
x = logspace(-2,-1.5,5);y = 1.5*x.^1;
plot(x,y,'k--','linewidth',2)
x = logspace(1.3,2,5);y = 0.7*x.^(1/2);
plot(x,y,'k--','linewidth',2)

xlabel('\epsilonN');
ylabel('R_{N}')
set(gca,'xscale','log');set(gca,'yscale','log');
box on
text(1.2e-2,3.5e-2,'N','FontSize',20)
text(40,3.5,'N^{1/2}','FontSize',20)

savename = sprintf('../../results/randcopoly-results/random-simulation/figure7B.eps');
saveas(gcf,savename,'epsc')

