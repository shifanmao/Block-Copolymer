function [KS_SIM,SINV_SIM,D2S_SIM,S,CHI]=plotsim5(EPS,LAM,vfac,PLOTON,SIMNUM,CHEMNUM)
% Plots density-density correlation of sim. and theory
% INPUTS::
%   EPS, number of Kuhn steps per monomer
%   LAM, degree of chemical correlation

% simulation folder
% folder = '../../results/randcopoly-results/scalcbatch-12-15-15';
folder = '../../results/randcopoly-results/scalcbatch-05-19-16';
addpath('../functions/')
addpath('../misc/')

% load simulation parameters
simparam=load([folder,'/chivals']);

if nargin == 4
    % find corresponding simulation at EPS and LAM
    ind = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
    % finding simulation index(indices)
    chemparam=load([folder,'/chemind']);
    SIMNUM = ind;
    CHEMNUM=chemparam(chemparam(:,1)==ind,2);
end

% simulation constants
FA=0.5;  % fraction of A blocks
M=8;        % number of blocks
EPS=EPS*vfac;
G=5*vfac;  % number of discrete monomers
NM=G*EPS;  % number of Kuhn steps per monomer

% range of chi params.
chiind=11;
plotind=chiind(1:end);

% Flory-Huggins parameter
chiv = load(sprintf([folder,'/sdata-%d-%d/Sdata/chilist'],SIMNUM,CHEMNUM));
chiv = chiv/G;

% Find spinodal CHI and critical k
[KS_MF,sval,d2gam2]=kmaxwlc(M,NM,FA,LAM);
chis=0.5*sval;

R2=-0.5+0.5*exp(-2*NM)+NM;
k=logspace(log10(0.01),log10(100),100)./sqrt(R2); % wavevectors

% Plot: density-density correlation
if PLOTON==1
    f1=figure;hold;set(gca,'fontsize',30)
    set(f1,'position',[0,0,800,600])
    cnt=find(chiind==plotind(1));
    for ii = plotind
        CHI = chiv(ii);
        col = cnt/length(chiind);

        if CHI<chis*EPS*0.85
            % Plot analytical theory (note conversion of monomer volume v, factor of EPS)
            val = s2invwlc(M,NM,FA,LAM,k);
            plot(k.*sqrt(R2),1./(-2*CHI+EPS*val),'--','color',[1-col 0 col],'linewidth',3)
        end

        % Plot simulation results
        filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
            SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
        S = load(filename);
        plot(S(:,1),S(:,2),'.-','color',[1-col 0 col],'linewidth',3,'markersize',20)
        cnt = cnt+1;
    end

    xlabel('R_Mq');
    ylabel('<\psi(q)\psi(-q)>')
    xlim([.6,12]);ylim([1e-1,2e2]);
    set(gca,'Ytick',[1e-1 1e0 1e1 1e2])
    set(gca,'Xtick',[1,10])
    set(gca,'xscale','log');set(gca,'yscale','log')
    box on
end

chiv=chiv(chiind);
% Find peak of structure factors
KS_SIM = zeros(length(chiv),1);
SINV_MF = zeros(length(chiv),1);
SINV_SIM = zeros(length(chiv),1);
D2S_SIM = zeros(length(chiv),1);
D2S_MF = zeros(length(chiv),1);
ERR = zeros(length(chiv),3);

for ii = 1:length(chiv)
% for ii = plotind
    CHI = chiv(ii);
    
    % Find peak position
    SINV_MF(ii)=-2*CHI+EPS*sval;
    D2S_MF(ii) = -1/(sval^2*R2*EPS)*d2gam2;

    % Plot simulation results
    filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
        SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
    S = load(filename);
    
    IND = find(S(:,2)==max(S(:,2)));IND = IND(1);
    KS_SIM(ii) = S(IND,1);
    SINV_SIM(ii) = 1./S(IND,2);
    
    NUMFIT = 4;
    if IND>NUMFIT  % central differences
        Kfit = S(IND-NUMFIT:IND+NUMFIT,1);
        Sfit = S(IND-NUMFIT:IND+NUMFIT,2);
    
        % local fit to Lorentzian (three parameter fit)
        options = optimset('Display','off',...
            'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e10,'MaxIter',1e10);
        x0 = [KS_SIM(ii),SINV_SIM(ii),1];   % x = [ks_sim,sinv_sim,d2s]
        lb = [0.1,0.01,0.1]; ub = [5,10,10];
        fun = @(x,Kfit) x(2) + (1/2)*(1/G)*x(3)*(Kfit-x(1)).^2;
        [xhat,~,resid,~,~,~,J] = lsqcurvefit(fun,x0,Kfit,1./Sfit,lb,ub,options);
        [~,sehat] = nlparcinew(xhat,resid,'jacobian',J);      % %estimate error from the fit
        x=xhat;se=sehat;

        %[x,r,J,cov,mse] = nlinfit(Kfit,1./Sfit,fun,x0,options);
        %[ci,se] = nlparcinew(x,r,'covar',cov);
    elseif IND<=NUMFIT  % forward differences
        Kfit = S(IND:IND+NUMFIT,1);
        Sfit = S(IND:IND+NUMFIT,2);
        
        % local fit to Lorentzian (two parameter fit)
        options = optimset('Display','off',...
            'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e10,'MaxIter',1e10);
        x0 = [SINV_SIM(ii),1];   % x = [sinv_sim,d2s]
        lb = [0.01,0.1]; ub = [10,10];
        fun = @(x,Kfit) x(1) + (1/2)*(1/G)*x(2)*(Kfit).^2;
        [xhat,~,resid,~,~,~,J] = lsqcurvefit(fun,x0,Kfit,1./Sfit,lb,ub,options);
        [~,sehat] = nlparcinew(xhat,resid,'jacobian',J);      % %estimate error from the fit
        x=[0,xhat];se=[0,sehat'];
    else
        error('peak location has large variability')
    end
    
    KS_SIM(ii) = x(1);
    SINV_SIM(ii) = x(2);
    D2S_SIM(ii) = x(3);
    ERR(ii,1:3) = se(1:3);

%     savename = sprintf('../../results/randcopoly-results/structure-figures/sfig-eps%.2f-lam%.2f.eps',EPS,LAM);
%     saveas(gcf,savename,'epsc')
end