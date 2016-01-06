% Figure 1: density correlations
clear;
CHAIN=2; % Gaussian chain if =1; WLC if =2

% Parameters
NM=0.1;  % number of Kuhn steps per monomer
N=100;  % number of monomers

FA=0.5;  % fraction of A monomers
LAM=-.75; % chemical correlation

ORDEig=20;  % number of eigvenvalues
NumLayer=500;  % number of residual layers

% find spinodal
[CHIS,ks,d2gam2]=spinodal(CHAIN,N,NM,LAM,FA,ORDEig,NumLayer);

% define wavevector
if CHAIN==1
    R2 = NM;
elseif CHAIN==2
    R2 = NM-(1/2)*(1-exp(-2*NM));
elseif CHAIN==3
    R2 = NM^2;
else
    error('Chain Definition error :: Gaussian Chain or Worm-like Chain or Rigid Rod')
end
k=logspace(-2,2,200)/sqrt(R2);

%% make plot
f1=figure;hold;set(gca,'fontsize',18);leg=[];
set(f1,'position',[0,0,800,600])
CHIV = fliplr([0.0:0.2:0.8]);
cnt=1;
col=flipud(copper(length(CHIV)));
col(:,2)=0;  % red to blue
for CHI = CHIV*CHIS/NM;
    CHI
    leg = [leg;sprintf('\\chi=%.1f \\chi_s',CHI/CHIS)];
    
    % calculate Gamma2
    G=gamma2(CHAIN,N,NM,LAM,FA,k,CHI,ORDEig,NumLayer);
    plot(k*sqrt(R2),1./(G),'color',col(cnt,:),'linewidth',2);
    cnt = cnt+1;
end

xlabel('kR_M');
ylabel('Density correlation $<\tilde{\psi}(\vec{k}) \tilde{\psi}(-\vec{k})>$','Interpreter','LaTex')
set(gca,'xscale','log');
if LAM<0
    legend(leg,'location','northwest')
else
    legend(leg,'location','northeast')
end
set(gca,'xscale','log')
set(gca,'yscale','log')