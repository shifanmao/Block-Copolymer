clear;close all;
addpath('functions')
addpath('chainstats')
addpath('misc')
addpath('chainstats/eigcalc')
addpath('chainstats/integrals')

N=1e1;  % number of statistical steps of total chain
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition
C=1e3;  % dimensionless excluded volume parameter in the Gaussian limit
        % In the Gaussian chain limit, Nbar = C^2

FIG1=0;     % Figure 1: make a mean-field phase diagram
FIG2=0;     % Figure 2: make a phase diagram with density fluctuations
FIG3=0;     % Figure 3: mean-field spinodal and critical wavelength at FA=0.5
FIG4=0;     % Figure 4: renormalized spinodal at FA=0.5
FIG5=1;     % Figure 5: density-density correlations
FIG6=0;     % Figure 6: vertex functions

% Figure 1: make a mean-field phase diagram
if (FIG1)
    plotphase(N,FAV);
end

% Figure 2: make a phase diagram with density fluctuations
if (FIG2)
    plotphaseRG(N,C,FAV);
end

% Figure 3: mean-field spinodal and critical wavelength at FA=0.5
if (FIG3)
    NV=logspace(-1,3,20);  % number of statistical steps of total chain
    [chis,ks,d2gam2]=spinodal(NV,0.5);
    figure;semilogx(NV,chis.*NV);xlabel('N');ylabel('\chiN')
    figure;loglog(NV,1./ks);xlabel('N');ylabel('1/q^*')
end

% Figure 4: renormalized spinodal at FA=0.5
if (FIG4)
    CV=logspace(1,4,21);
    [chit,phase]=spinodalRG(N,CV,0.5);
    chit=reshape(chit,length(CV),1);

    figure;hold;set(gca,'fontsize',20)
    col='b';
    plot(CV.^2,ones(length(CV),1)*spinodal(N,0.5)*N,'--','linewidth',2,'color',col)
    plot(CV.^2,chit*N,'s','MarkerSize',8,'MarkerFaceColor',col,'MarkerEdgeColor',col);

    % %Empirical solutions
    plot(CV.^2,ones(length(CV),1)*10.495,'--','linewidth',2,'color','k')
    set(gca,'xscale','log');box on
    xlabel('C^2');ylabel('\chiN');title(['N=',num2str(N)])
    legend('MF theory','Renormalized ODT','Fit')
end

% Figure 5: density-density correlations
if (FIG5)
    densityRG(N,C,0.5);
end

% Figure 6: vertex functions
if (FIG6)
    NQ=1;  % number of wavevector sets in calculating GAM4
    k = logspace(-1,2,100);
    G=gamma2(N,FA,k,0);
    figure;plot(k,1./G);
    [gam3,gam4]=calcgamma(N,FAV,NQ);
    figure;plot(FAV,-gam3*N,'k-','linewidth',2);xlim([0.2,0.5]);
    xlabel('f_A');ylabel('-N\Gamma_3(q^*)')
    figure;plot(FAV,gam4*N,'k-','linewidth',2);xlim([0.3,0.5]);
    xlabel('f_A');ylabel('N\Gamma_4(q^*)')
end

% Save to images
% saveas(gca,'../results/diblockcopoly-results/figures/phaseRGC1e2_flexible.eps','epsc')