function [chis,chit]=densityRG(N,Nbar,FA)
% Plot density-density correlations during phase transition
% according to Mean-field theory and FH theory
% Usage :: [chis,chit]=densityRG(N,Nbar,FA)
% Inputs ::
%   N, number of Kuhn steps of whole chain
%   Nbar, invariant degree of polymerization
%   FA, fraction of A monomers
%   CHIV, range of CHI values

% find mean-field spinodal
[chis,ks,d2gamma2]=spinodal(N,FA);

% PLOT1: DENSITY CORRELATION  %%
% wavevectors
kmin=0.5;kmax=3.5;
k=power(linspace(kmin,kmax,200),1)/sqrt(N/6);

% Flory-Huggins parameter
CHIV=linspace(0,1,5);

figure;hold;set(gca,'fontsize',20)
for ii = 1:length(CHIV)
    CHI=chis*CHIV(ii);
    if length(CHIV)>1
        col = (ii-1)/(length(CHIV)-1);
    else
        col = 0;
    end

    if CHI<=0.9*chis
        % plot mean-field results
        Smf=densitymf(N,FA,k,CHI);
        plot(k*sqrt(N/6),Smf./N,'color',[col 0 1-col],'linestyle','--','linewidth',2);

        % plot RG results
        Sfh=densityfh(N,Nbar,FA,k,ks,CHI,d2gamma2);
        plot(k*sqrt(N/6),Sfh./N,'color',[col 0 1-col],'linestyle','-','linewidth',2);
    end
end
xlim([kmin,kmax]);box on
xlabel('qR');ylabel('<\psi^2>/N')

% PLOT2: CRITICAL MODE vs CHI
% Flory-Huggins parameter
CHIV=linspace(0,4,50);

Smf = zeros(length(CHIV),1);
Sfh = zeros(length(CHIV),1);
for ii = 1:length(CHIV)
    CHI = CHIV(ii)*chis;
    Smf(ii)=densitymf(N,FA,ks,CHI);
    Sfh(ii)=densityfh(N,Nbar,FA,ks,ks,CHI,d2gamma2);
end

figure;hold;set(gca,'fontsize',20)
plot(CHIV*chis*N,N./Smf,'k--','linewidth',2);
plot(CHIV*chis*N,N./Sfh,'k-','linewidth',2);
% plot(CHIV*N,N./Swlc,'k-','linewidth',2);
xlim([1,17]);ylim([0,20]);box on
xlabel('\chi N');ylabel('N/<\psi^2(q^*)>')

% find renormalized spinodal
chit=spinodalRG(N,Nbar,FA);
plot(chit*N,N./densityfh(N,Nbar,FA,ks,ks,chit,d2gamma2),'ks',...
    'MarkerSize',6,'MarkerFaceColor','k');
end

function Smf=densitymf(N,FA,k,CHI)
% Mean-field theory (Leibler)
Gmf=gamma2(N,FA,k,CHI);
Smf=1./Gmf;
end
function Sfh=densityfh(N,Nbar,FA,k,ks,CHI,d2gamma2)
% RG with no q dependence (F-H)
% self-consistent equations:
%  unknowns x(1)=r, x(2)=alpha
%  r = N*gamma2+pref*N*gamma4
%  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
%       , where pref = N*Qs/(4*pi*power(Nbar*r*alpha,1/2))

% Calculate a few constants gamma2
gam2 = gamma2(N,FA,ks,CHI);

% gamma4 at angle phi=pi (assume no q dependence)
NQ=1;
[~,gam4]=calcgamma(N,FA,NQ);
gam4=real(gam4(end,1));

% solve self-consistent equations
alpha = (1/2)*d2gamma2;
pref = (N*ks^2)/(4*pi*power(Nbar*alpha,1/2));
root = roots([1,0,-N*gam2,-pref*N*gam4]);
r = power(root(root>0),2);

Gfh = r/N + alpha*(k-ks).^2;
Sfh=1./(Gfh);
end

% Reserved for future work
% %% RG with no q dependence (F-H)
% CHAIN=1;
% ORDEig=4;NumLayer=500;
% % self-consistent equations:
% %  unknowns x(1)=r, x(2)=alpha
% %  r = N*gamma2+pref*N*gamma4
% %  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
% %       , where pref = N*Qs/(4*pi*power(Nbar*r*alpha,1/2))
% 
% % Calculate a few constants
% % gamma2
% gam2 = gamma2(CHAIN,N,NM,FA,Qs,0,ORDEig,NumLayer);
% % second-order derivative of gamma2 near Qs
% dQ=Qs*1e-3;
% d2gamma2=(gamma2(CHAIN,N,NM,FA,Qs+dQ,0,ORDEig,NumLayer)-...
%             2*gamma2(CHAIN,N,NM,FA,Qs,0,ORDEig,NumLayer)+...
%             gamma2(CHAIN,N,NM,FA,Qs-dQ,0,ORDEig,NumLayer))/(dQ^2);
% % gamma4 at angle phi=pi (assume no q dependence)
% Q1=[1,0,0]';Q2=rotz(pi)*Q1;Q3=-Q2;Q4=-Q1;
% gam4 = gamma4(Qs,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer);
% d2gamma4 = (gamma4(Qs+dQ,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer)-...
%             2*gamma4(Qs,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer)+...
%             gamma4(Qs-dQ,Q1,Q2,Q3,Q4,CHAIN,N,NM,FA,ORDEig,ORDEig,NumLayer))/(dQ^2);
%         
% % % initial guesses, lower bounds and upper bounds of solutions
% % unknowns: x = [r,alpha]
% x0 = [1,1];
% lb=[0,0];
% ub=[1e5,1e5];
% options = optimset('Display','off',...
%     'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);
% F = @(x) [x(1)-N*gam2-(N*Qs^2)/(4*pi*power(Nbar*x(1)*x(2),1/2))*N*gam4,...
%           2*x(2)-d2gamma2-(N*Qs^2)/(4*pi*power(Nbar*x(1)*x(2),1/2))*d2gamma4];
% [x,~] = lsqnonlin(F,x0,lb,ub,options);
% r = x(1); alpha = x(2);
% Gfh = r/N + alpha*(Q-Qs).^2;
% chisfh = (r/N)/2*N

% % RG with q dependence (Fredrickson-Barret)
% % self-consistent equations:
% %  unknowns x(1)=r, x(2)=Qstar, x(3)=alpha
% %  r = N*gamma2+pref*N*gamma4
% %  0 = d(gamma2)/dQ + pref*d(gamma4)/dQ
% %  2*alpha = d^2(gamma2)/dQ^2 + pref*d^2(gamma4)/dQ^2
% %
% %  where pref = N*Qstar^2/(4*pi*power(Nbar*r*alpha,1/2))
% r = 0; alpha = 0; Qstar = Qs;
%     
% % initial guesses, lower bounds and upper bounds of solutions
% x0 = [1,Qs,1];
% lb=[0,0,0];
% ub=[1e5,1e5,1e5];
% options = optimset('Display','off',...
%     'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);
% F = @(x) [x(1)-N*gamma2-(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*N*gamma4,...
%           dgamma2+(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*dgamma4,...
%           x(3)-d2gamma2-(N*x(2)^2)/(4*pi*power(Nbar*x(1)*x(3),1/2))*d2gamma4];
% [x,~] = lsqnonlin(F,x0,lb,ub,options);
% G = r/N + alpha*(Q-Qstar).^2;