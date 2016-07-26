function [F1,F3,F6] = testgam(NV)
%clear;
addpath('../functions/')
addpath('../chainstats/')
addpath('../chainstats/eigcalc/')
addpath('../chainstats/integrals/')
addpath('../misc/')

%NV=logspace(-1,4,21);
%NV = 1e1;
FAV=linspace(0.1,0.5,41);  % range of A monomer chemical composition

%alpha=power(d2gam2/2.*NV./r2(NV),1/2);

NQ=4;
[gam3,gam4]=calcgamma(NV,FAV,NQ);

% % Plot Gammas
% figure;plot(FAV,-gam3*NV,'linewidth',2);xlim([.2,.5])
% set(gca,'fontsize',20);xlabel('f_A');ylabel('-N\Gamma^{(3)}')
% 
% figure;plot(FAV,gam4*NV,'linewidth',2);xlim([.3,.5])
% legend('\theta=0','\theta=\pi/3','\theta=pi/2','(1,2)')
% set(gca,'fontsize',20);xlabel('f_A');ylabel('N\Gamma^{(4)}')

% Plot cofficients
alpha1 = zeros(length(gam3),1);
beta1 = (1/4)*gam4(:,1);

alpha3 = (2/3/sqrt(3))*gam3;
beta3 = (1/12)*(gam4(:,1)+4*gam4(:,2));

alpha6 = (4/3/sqrt(6))*gam3;
beta6 = (1/24)*(gam4(:,1)+8*gam4(:,2)+2*gam4(:,3)+4*gam4(:,4));

figure;plot(FAV,-[alpha1,alpha3,alpha6]*NV);
legend('LAM','HEX','BCC')
set(gca,'fontsize',20);xlabel('f_A');ylabel('\alpha')

figure;plot(FAV,[beta1,beta3,beta6]*NV);
legend('LAM','HEX','BCC')
set(gca,'fontsize',20);xlabel('f_A');ylabel('\beta')

figure;hold;
CHI = 1/NV;
F1 = zeros(length(FAV),1);
F3 = zeros(length(FAV),1);
F6 = zeros(length(FAV),1);
for ii = 1:length(FAV)
    %[chis,ks,d2gam2]=spinodal(NV,FAV(ii));
    
    F1(ii) = lamellar(CHI,gam4(ii,1:4));
    F3(ii) = hexagonal(CHI,gam3(ii),gam4(ii,1:4));
    F6(ii) = bcc(CHI,gam3(ii),gam4(ii,1:4));
end
plot(FAV,F1*NV,'k--','linewidth',2)
plot(FAV,F3*NV,'b--','linewidth',2)
plot(FAV,F6*NV,'r--','linewidth',2)

CHI = 10/NV;
F1 = zeros(length(FAV),1);
F3 = zeros(length(FAV),1);
F6 = zeros(length(FAV),1);
for ii = 1:length(FAV)
    %[chis,ks,d2gam2]=spinodal(NV,FAV(ii));
    
    F1(ii) = lamellar(CHI,gam4(ii,1:4));
    F3(ii) = hexagonal(CHI,gam3(ii),gam4(ii,1:4));
    F6(ii) = bcc(CHI,gam3(ii),gam4(ii,1:4));
end
plot(FAV,F1*NV,'k-','linewidth',2)
plot(FAV,F3*NV,'b-','linewidth',2)
plot(FAV,F6*NV,'r-','linewidth',2)
legend('LAM','HEX','BCC')
set(gca,'yscale','log')
set(gca,'fontsize',20);xlabel('f_A');ylabel('\betaF');box on

end

function F1=lamellar(CHI,gam4)
    %%% LAM phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B1=0;
    C1=power(2*pi,-9)*(1/4)*gam4(1);

    % free energy %
    FE=@(psi) A*power(psi,2)+B1*power(psi,3)+C1*power(psi,4);
    x1=fminbnd(FE,1e-10,1e5);
    F1=FE(x1);
end
function F3=hexagonal(CHI,gam3,gam4)
    %%% HEX phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B3=power(2*pi,-6)*(2/3/sqrt(3))*gam3;
    C3=power(2*pi,-9)*(1/12)*(gam4(1)+4*gam4(2));

    % free energy %
    FE=@(psi) A*power(psi,2)+B3*power(psi,3)+C3*power(psi,4);
    x3=fminbnd(FE,1e-10,1e5);
    F3=FE(x3);
end
function F6=bcc(CHI,gam3,gam4)
    %%% BCC phase %%%
    % coefficients %
    A=-power(2*pi,-3)*2*CHI;
    B6=power(2*pi,-6)*(4/3/sqrt(6))*gam3;
    C6=power(2*pi,-9)*(1/24)*(gam4(1)+8*gam4(2)+2*gam4(3)+4*gam4(4));

    % free energy %
    FE=@(psi) A*power(psi,2)+B6*power(psi,3)+C6*power(psi,4);
    x6=fminbnd(FE,1e-10,1e5);
    F6=FE(x6);
end