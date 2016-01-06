% fhphase.m :: This function predicts diblock copolymer phase diagram
% with Fredrickson-Helfand density fluctuation correction
clear all;

readin=1;

NMV=1:6;
CHAIN=2;            % Which type of chain? 1=Gaussian, 2=WLC, 3=Rigid Rod
NMV=0;

for ii=1:length(NMV)
NM=10^NMV(ii)
NM=5;

% parameters for Worm-like chain propagator calculations
ORDEig=4;          % Number of eigenvalues
ORDL=5;             % Number of spherical harmonics
NumLayer=500;     % Number of residual layers

% parameters
FAV=linspace(0.1,0.5,21);               % fraction of A-type monomers
NFA=length(FAV);
NQ=1;

% results to return
chiS=zeros(NFA,1);      % spinodal (uncorrected)
chit=zeros(NFA,1);      % spinodal (corrected)
chi1=zeros(NFA,1);      % spinodal (corrected)
chi3=zeros(NFA,1);      % spinodal (corrected)
chi6=zeros(NFA,1);      % spinodal (corrected)
chi13=zeros(NFA,1);      % spinodal (corrected)
chi36=zeros(NFA,1);      % spinodal (corrected)
ks=zeros(NFA,1);

% parameters
c=zeros(NFA,1);
d=zeros(NFA,1);
d2c=zeros(NFA,1);
    
% options = optimset('Display','off',...
%     'TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',1e14,'MaxIter',1e14);
options = optimset('Display','off',...
    'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);

if (readin)
    filename=sprintf('data/gam_C%dNM%.2f.mat',CHAIN,NM);
    load(filename);
    gamma3 = gam3;
    gamma4 = gam4;
    NFA=length(FAV);
else
    gamma3 = zeros(length(FAV),1);
    gamma4 = zeros(length(FAV),1);
end

for ii=1:NFA
    FA=FAV(ii);
    
    % initial guesses, lower bounds and upper bounds of solutions
    x0 = [-1e2,1,1,1];
    lb=[-1e7,0,0,0];
    ub=[0,1e5,1e5,1e5];

    x02 = [-1e2,1,1,1,1,1,1];
    lb2=[-1e5,0,0,0,0,0,0];
    ub2=[0,1e5,1e5,1e5,1e5,1e5,1e5];
    
    if (readin)
        FA        
        chiS(ii)=chis(ii);
        gam3 = real(gamma3(ii));
        gam4 = real(gamma4(ii,1));
    else
    
        % calculate vertex functions
        if CHAIN==1
            [chiS(ii),ks(ii),gam3,gam4]=calcgamma(CHAIN,NM,FA,ORDEig,ORDL,NumLayer,NQ);
        elseif CHAIN==2
            [chiS(ii),ks(ii),gam3,gam4]=calcgamma(CHAIN,NM,FA,ORDEig,ORDL,NumLayer,NQ);
        end
        gam3=real(gam3);
        gam4=real(gam4);
        
        gamma3(ii)=gam3;
        gamma4(ii)=gam4;
    end
    
    % calculate constant (estimate local second-order derivative)
    xs=ks(ii)^2*(1/6)*NM;
    dx=xs*1e-3;
    dks=dx/(2*ks(ii)*NM/6);
    d2c(ii)=NM*(gamma2(CHAIN,NM,FA,ks(ii)+dks,0,20,NumLayer)-...
                2*gamma2(CHAIN,NM,FA,ks(ii),0,20,NumLayer)+...
                gamma2(CHAIN,NM,FA,ks(ii)-dks,0,20,NumLayer))/(dx^2);
	c(ii)=power(1/3*xs*d2c(ii),1/2);
    
    % parameters
    miu = NM*gam3/power(c(ii),3);
    lam = NM*gam4(1)/power(c(ii),4);
    d(ii) = (3*xs)/(2*pi);

    %%%%% LAM phase %%%%%
    n = 1;
    if n==1
        theta=0*miu;
        eta=-lam/2;
    elseif n==3
        theta=-miu;
        eta=-lam/2;
    elseif n==6
        theta=-2*miu;
        eta=3*lam/2;
    end
    
    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0, x(3)=r, x(4)=a
    % Equation 1: r0 = tau+d*lam*power(r0*NM,-1/2)
    % Equation 2: r = tau+d*lam*power(r*NM,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(NM,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*NM,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*NM,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(NM,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi1(ii)=chiS(ii)-c(ii)^2*tau/(2*NM);
    
    %%%%% HEX phase %%%%%
    n = 3;
    if n==1
        theta=0*miu;
        eta=-lam/2;
    elseif n==3
        theta=-miu;
        eta=-lam/2;
    elseif n==6
        theta=-2*miu;
        eta=3*lam/2;
    end
    
    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0, x(3)=r, x(4)=a
    % Equation 1: r0 = tau+d*lam*power(r0*NM,-1/2)
    % Equation 2: r = tau+d*lam*power(r*NM,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(NM,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*NM,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*NM,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(NM,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi3(ii)=chiS(ii)-c(ii)^2*tau/(2*NM);
    
    %%%%% BCC phase %%%%%
    n = 6;
    if n==1
        theta=0*miu;
        eta=-lam/2;
    elseif n==3
        theta=-miu;
        eta=-lam/2;
    elseif n==6
        theta=-2*miu;
        eta=3*lam/2;
    end
    
    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0, x(3)=r, x(4)=a
    % Equation 1: r0 = tau+d*lam*power(r0*NM,-1/2)
    % Equation 2: r = tau+d*lam*power(r*NM,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(NM,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*NM,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*NM,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(NM,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi6(ii)=chiS(ii)-c(ii)^2*tau/(2*NM);
    chit(ii)=min([chi1(ii),chi3(ii),chi6(ii)]);

    %%%%% LAM/HEX phase %%%%%

    if (chit(ii)==chi6(ii) || chit(ii)==chi3(ii))
        theta_1=0*miu;
        eta_1=-lam/2;
        theta_3=-miu;
        eta_3=-lam/2;

        % find corrected spinodal
        % use Hartree approximation from Fredrickson-Helfand theory
        % to evalute free energy of diblock copolymers

        % unknowns x(1)=tau, x(2)=r0_1(LAM), x(3)=r_1(LAM), x(4)=a_1(LAM)
        % unknowns           x(5)=r0_3(HEX), x(6)=r_3(HEX), x(7)=a_3(HEX)
        % Equation 1: r0_1 = tau+d*lam*power(r0_1*NM,-1/2) (LAM)
        % Equation 2: r0_3 = tau+d*lam*power(r0_3*NM,-1/2) (HEX)
        % Equation 3: r_1 = tau+d*lam*power(r_1*NM,-1/2)+1*lam*a_1^2 (LAM)
        % Equation 4: r_3 = tau+d*lam*power(r_3*NM,-1/2)+3*lam*a_3^2 (HEX)
        % Equation 5: eta_1*a_1^2-theta_1*a+r_1 = 0
        % Equation 6: eta_3*a_3^2-theta_3*a+r_3 = 0
        % Equation 7: phi_1 = phi_3 = (1/2/lam)*(r_1^2-r0_1^2)+...
        %               d*power(NM,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4 - ...
        %                             (1/2/lam)*(r_1^3-r0_1^3)+...
        %               d*power(NM,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4

        F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*NM,-1/2),...
                  x(5)-x(1)-d(ii)*lam*power(x(5)*NM,-1/2),...
                  x(3)-x(1)-d(ii)*lam*power(x(3)*NM,-1/2)-1*lam*x(4)^2,...
                  x(6)-x(1)-d(ii)*lam*power(x(6)*NM,-1/2)-3*lam*x(7)^2,...
                  eta_1*x(4)^2-theta_1*x(4)+x(3),...
                  eta_3*x(7)^2-theta_3*x(7)+x(6),...
                  ((1/2/lam)*(x(3)^2-x(2)^2)+d(ii)*power(NM,-1/2)*(sqrt(x(3))-sqrt(x(2)))-...
                    2*1/3*theta_1*x(4)^3+(1/2)*1*eta_1*x(4)^4)-...
                  ((1/2/lam)*(x(6)^2-x(5)^2)+d(ii)*power(NM,-1/2)*(sqrt(x(6))-sqrt(x(5)))-...
                    2*3/3*theta_3*x(7)^3+(1/2)*3*eta_3*x(7)^4)];
        [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
        tau=x(1);
        chi13(ii)=chiS(ii)-c(ii)^2*tau/(2*NM);
    else
        chi13(ii)=0;
    end
    
    %%%%% HEX/BCC phase %%%%%
    if (chit(ii)==chi6(ii))
        theta_3=-miu;
        eta_3=-lam/2;
        theta_6=-2*miu;
        eta_6=3*lam/2;

        % find corrected spinodal
        % use Hartree approximation from Fredrickson-Helfand theory
        % to evalute free energy of diblock copolymers

        % unknowns x(1)=tau, x(2)=r0_1(HEX), x(3)=r_1(HEX), x(4)=a_1(HEX)
        % unknowns           x(5)=r0_3(BCC), x(6)=r_3(BCC), x(7)=a_3(BCC)
        % Equation 1: r0_1 = tau+d*lam*power(r0_1*NM,-1/2) (HEX)
        % Equation 2: r0_3 = tau+d*lam*power(r0_3*NM,-1/2) (BCC)
        % Equation 3: r_1 = tau+d*lam*power(r_1*NM,-1/2)+1*lam*a_1^2 (HEX)
        % Equation 4: r_3 = tau+d*lam*power(r_3*NM,-1/2)+3*lam*a_3^2 (BCC)
        % Equation 5: eta_1*a_1^2-theta_1*a+r_1 = 0
        % Equation 6: eta_3*a_3^2-theta_3*a+r_3 = 0
        % Equation 7: phi_1 = phi_3 = (1/2/lam)*(r_1^2-r0_1^2)+...
        %               d*power(NM,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4 - ...
        %                             (1/2/lam)*(r_1^3-r0_1^3)+...
        %               d*power(NM,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4

        F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*NM,-1/2),...
                  x(5)-x(1)-d(ii)*lam*power(x(5)*NM,-1/2),...
                  x(3)-x(1)-d(ii)*lam*power(x(3)*NM,-1/2)-3*lam*x(4)^2,...
                  x(6)-x(1)-d(ii)*lam*power(x(6)*NM,-1/2)-6*lam*x(7)^2,...
                  eta_3*x(4)^2-theta_3*x(4)+x(3),...
                  eta_6*x(7)^2-theta_6*x(7)+x(6),...
                  ((1/2/lam)*(x(3)^2-x(2)^2)+d(ii)*power(NM,-1/2)*(sqrt(x(3))-sqrt(x(2)))-...
                    2*3/3*theta_3*x(4)^3+(1/2)*3*eta_3*x(4)^4)-...
                  ((1/2/lam)*(x(6)^2-x(5)^2)+d(ii)*power(NM,-1/2)*(sqrt(x(6))-sqrt(x(5)))-...
                    2*6/3*theta_6*x(7)^3+(1/2)*6*eta_6*x(7)^4)];
        [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
        tau=x(1);
        chi36(ii)=chiS(ii)-c(ii)^2*tau/(2*NM);
    else
        chi36(ii)=0;
    end
end

ind13=find(chit<chi1);
ind36=find(chit==chi6);
% ind13=1:NFA;ind36=1:NFA;
figure;hold;set(gca,'fontsize',18)
plot(FAV,chiS*NM,'k--','linewidth',1.5)
plot(FAV,chit*NM,'k','linewidth',1.5)
plot(1-FAV,chiS*NM,'k--','linewidth',1.5)
plot(1-FAV,chit*NM,'k','linewidth',1.5)
plot(FAV(ind13),chi13(ind13)*NM,'r','linewidth',1.2)
plot(FAV(ind36),chi36(ind36)*NM,'b','linewidth',1.2)
plot(1-FAV(ind13),chi13(ind13)*NM,'r','linewidth',1.2)
plot(1-FAV(ind36),chi36(ind36)*NM,'b','linewidth',1.2)
xlabel('f');ylabel('\chi N')
xlim([FAV(1),1-FAV(1)])

figure;hold;set(gca,'fontsize',18)
plot(FAV,chiS*NM,'k--','linewidth',1.5)
plot(FAV,chi1*NM,'k-','linewidth',1.5)
plot(FAV,chi3*NM,'r.','linewidth',1.5)
plot(FAV,chi6*NM,'b-','linewidth',1.5)
xlabel('f');ylabel('\chi N')

% save data
if (NM>=10 && mod(NM,10)==0)
    filename=sprintf('data/phase_C%dNM1e%d.mat',CHAIN,log10(NM))
    save(filename,'FAV','ks','chiS','chit','chi1','chi3','chi6','chi13','chi36','gamma3','gamma4','c','d','d2c')
else
    filename=sprintf('data/phase_C%dNM%.2f.mat',CHAIN,NM)
    save(filename,'FAV','ks','chiS','chit','chi1','chi3','chi6','chi13','chi36','gamma3','gamma4','c','d','d2c')
end
end