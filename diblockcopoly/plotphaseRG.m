function plotphaseRG(N,Nbar,FAV)
% fhphase.m :: This function predicts diblock copolymer phase diagram
% with Fredrickson-Helfand density fluctuation correction
% Usage :: [chis,chi13,chi36,chi12,chi23,chi26]=plotphaseRG(N,Nbar,FAV)
% Inputs ::
%    FAV, fraction of A-type monomers

% results to return
NFA=length(FAV);
chis=zeros(NFA,1);      % spinodal (uncorrected)
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
gamma3 = zeros(NFA,1);
gamma4 = zeros(NFA,1);

for ii=1:NFA
    FA=FAV(ii);
    
    % calculate mean-field solution
    [chis(ii),ks(ii),~]=spinodal(N,FA);
        
    % calculate vertex functions
    NQ=1;  % assume to Q dependence
    [gam3,gam4]=calcgamma(N,FA,NQ);
    gam3=real(gam3);
    gam4=real(gam4);

    gamma3(ii)=gam3;
    gamma4(ii)=gam4;
    
    % calculate constant (estimate local second-order derivative)
    xs=ks(ii)^2*(1/6)*N;
    dx=xs*1e-3;
    dks=dx/(2*ks(ii)*N/6);
    d2c(ii)=N*(gamma2(N,FA,ks(ii)+dks,0)-2*gamma2(N,FA,ks(ii),0)+gamma2(N,FA,ks(ii)-dks,0))/(dx^2);
	c(ii)=power(1/3*xs*d2c(ii),1/2);
    
    % parameters
    miu = N*gam3/power(c(ii),3);
    lam = N*gam4(1)/power(c(ii),4);
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
    
    % Start solving self-consistent equations
    fprintf('Step 2: Calculating renormalized OOT phase diag. at FA=%.2f, N=%.2e\n',FA,N)
    
    % solver options
    options = optimset('Display','off',...
        'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);

    % initial guesses, lower bounds and upper bounds of solutions
    x0=[-1e2,1,1,1];
    lb=[-1e7,0,0,0];
    ub=[0,1e5,1e5,1e5];

    x02=[-1e2,1,1,1,1,1,1];
    lb2=[-1e5,0,0,0,0,0,0];
    ub2=[0,1e5,1e5,1e5,1e5,1e5,1e5];
    
    % find corrected spinodal
    % use Hartree approximation from Fredrickson-Helfand theory
    % to evalute free energy of diblock copolymers

    % unknowns x(1)=tau, x(2)=r0, x(3)=r, x(4)=a
    % Equation 1: r0 = tau+d*lam*power(r0*Nbar,-1/2)
    % Equation 2: r = tau+d*lam*power(r*Nbar,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(Nbar,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*Nbar,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*Nbar,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(Nbar,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi1(ii)=chis(ii)-c(ii)^2*tau/(2*N);
    
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
    % Equation 1: r0 = tau+d*lam*power(r0*Nbar,-1/2)
    % Equation 2: r = tau+d*lam*power(r*Nbar,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(Nbar,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*Nbar,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*Nbar,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(Nbar,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi3(ii)=chis(ii)-c(ii)^2*tau/(2*N);
    
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
    % Equation 1: r0 = tau+d*lam*power(r0*Nbar,-1/2)
    % Equation 2: r = tau+d*lam*power(r*Nbar,-1/2)+n*lam*a^2
    % Equation 3: eta*a^2-theta*a+r = 0
    % Equation 4: phi = (1/2/lam)*(r^2-r0^2)+...
    %               d*power(Nbar,-1/2)*(r^0.5-r0^0.5)...
    %               -2*n/3*theta*a^3+(1/2)*n*eta*a^4 = 0

    F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*Nbar,-1/2),...
              x(3)-x(1)-d(ii)*lam*power(x(3)*Nbar,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d(ii)*power(Nbar,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chi6(ii)=chis(ii)-c(ii)^2*tau/(2*N);
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
        % Equation 1: r0_1 = tau+d*lam*power(r0_1*Nbar,-1/2) (LAM)
        % Equation 2: r0_3 = tau+d*lam*power(r0_3*Nbar,-1/2) (HEX)
        % Equation 3: r_1 = tau+d*lam*power(r_1*Nbar,-1/2)+1*lam*a_1^2 (LAM)
        % Equation 4: r_3 = tau+d*lam*power(r_3*Nbar,-1/2)+3*lam*a_3^2 (HEX)
        % Equation 5: eta_1*a_1^2-theta_1*a+r_1 = 0
        % Equation 6: eta_3*a_3^2-theta_3*a+r_3 = 0
        % Equation 7: phi_1 = phi_3 = (1/2/lam)*(r_1^2-r0_1^2)+...
        %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4 - ...
        %                             (1/2/lam)*(r_1^3-r0_1^3)+...
        %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4

        F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*Nbar,-1/2),...
                  x(5)-x(1)-d(ii)*lam*power(x(5)*Nbar,-1/2),...
                  x(3)-x(1)-d(ii)*lam*power(x(3)*Nbar,-1/2)-1*lam*x(4)^2,...
                  x(6)-x(1)-d(ii)*lam*power(x(6)*Nbar,-1/2)-3*lam*x(7)^2,...
                  eta_1*x(4)^2-theta_1*x(4)+x(3),...
                  eta_3*x(7)^2-theta_3*x(7)+x(6),...
                  ((1/2/lam)*(x(3)^2-x(2)^2)+d(ii)*power(Nbar,-1/2)*(sqrt(x(3))-sqrt(x(2)))-...
                    2*1/3*theta_1*x(4)^3+(1/2)*1*eta_1*x(4)^4)-...
                  ((1/2/lam)*(x(6)^2-x(5)^2)+d(ii)*power(Nbar,-1/2)*(sqrt(x(6))-sqrt(x(5)))-...
                    2*3/3*theta_3*x(7)^3+(1/2)*3*eta_3*x(7)^4)];
        [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
        tau=x(1);
        chi13(ii)=chis(ii)-c(ii)^2*tau/(2*N);
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
        % Equation 1: r0_1 = tau+d*lam*power(r0_1*Nbar,-1/2) (HEX)
        % Equation 2: r0_3 = tau+d*lam*power(r0_3*Nbar,-1/2) (BCC)
        % Equation 3: r_1 = tau+d*lam*power(r_1*Nbar,-1/2)+1*lam*a_1^2 (HEX)
        % Equation 4: r_3 = tau+d*lam*power(r_3*Nbar,-1/2)+3*lam*a_3^2 (BCC)
        % Equation 5: eta_1*a_1^2-theta_1*a+r_1 = 0
        % Equation 6: eta_3*a_3^2-theta_3*a+r_3 = 0
        % Equation 7: phi_1 = phi_3 = (1/2/lam)*(r_1^2-r0_1^2)+...
        %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4 - ...
        %                             (1/2/lam)*(r_1^3-r0_1^3)+...
        %               d*power(Nbar,-1/2)*(r_1^0.5-r0_1^0.5)...
        %               -2*1/3*theta_1*a^3+(1/2)*1*eta_1*a^4

        F = @(x) [x(2)-x(1)-d(ii)*lam*power(x(2)*Nbar,-1/2),...
                  x(5)-x(1)-d(ii)*lam*power(x(5)*Nbar,-1/2),...
                  x(3)-x(1)-d(ii)*lam*power(x(3)*Nbar,-1/2)-3*lam*x(4)^2,...
                  x(6)-x(1)-d(ii)*lam*power(x(6)*Nbar,-1/2)-6*lam*x(7)^2,...
                  eta_3*x(4)^2-theta_3*x(4)+x(3),...
                  eta_6*x(7)^2-theta_6*x(7)+x(6),...
                  ((1/2/lam)*(x(3)^2-x(2)^2)+d(ii)*power(Nbar,-1/2)*(sqrt(x(3))-sqrt(x(2)))-...
                    2*3/3*theta_3*x(4)^3+(1/2)*3*eta_3*x(4)^4)-...
                  ((1/2/lam)*(x(6)^2-x(5)^2)+d(ii)*power(Nbar,-1/2)*(sqrt(x(6))-sqrt(x(5)))-...
                    2*6/3*theta_6*x(7)^3+(1/2)*6*eta_6*x(7)^4)];
        [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
        tau=x(1);
        chi36(ii)=chis(ii)-c(ii)^2*tau/(2*N);
    else
        chi36(ii)=0;
    end
end

ind13=find(chit<chi1);
ind36=find(chit==chi6);
% ind13=1:NFA;ind36=1:NFA;
figure;hold;set(gca,'fontsize',18)
plot(FAV,chis*N,'k--','linewidth',1.5)
plot(FAV,chit*N,'k','linewidth',1.5)
plot(1-FAV,chis*N,'k--','linewidth',1.5)
plot(1-FAV,chit*N,'k','linewidth',1.5)
plot(FAV(ind13),chi13(ind13)*N,'r','linewidth',1.2)
plot(FAV(ind36),chi36(ind36)*N,'b','linewidth',1.2)
plot(1-FAV(ind13),chi13(ind13)*N,'r','linewidth',1.2)
plot(1-FAV(ind36),chi36(ind36)*N,'b','linewidth',1.2)
xlabel('f');ylabel('\chi N')
xlim([FAV(1),1-FAV(1)])

figure;hold;set(gca,'fontsize',18)
plot(FAV,chis*N,'k--','linewidth',1.5)
plot(FAV,chi1*N,'k-','linewidth',1.5)
plot(FAV,chi3*N,'r.','linewidth',1.5)
plot(FAV,chi6*N,'b-','linewidth',1.5)
xlabel('f');ylabel('\chi N')

% % save data
% if (NM>=10 && mod(NM,10)==0)
%     filename=sprintf('data/phase_C%dNM1e%d.mat',CHAIN,log10(NM))
%     save(filename,'FAV','ks','chis','chit','chi1','chi3','chi6','chi13','chi36','gamma3','gamma4','c','d','d2c')
% else
%     filename=sprintf('data/phase_C%dNM%.2f.mat',CHAIN,NM)
%     save(filename,'FAV','ks','chis','chit','chi1','chi3','chi6','chi13','chi36','gamma3','gamma4','c','d','d2c')
% end