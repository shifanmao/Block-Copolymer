function [chit,phase,chi13,chi36]=plotphaseRG(N,Nbar,FAV)
% fhphase.m :: This function predicts diblock copolymer phase diagram
% with Fredrickson-Helfand density fluctuation correction
% Usage :: [chit,chi13,chi36,chi1,chi3,chi6]=plotphaseRG(N,Nbar,FAV)
% Inputs ::
%    FAV, fraction of A-type monomers

% results to return
NFA=length(FAV);
chis=zeros(NFA,1);      % spinodal (uncorrected)
chit=zeros(NFA,1);      % spinodal (corrected)
phase=zeros(NFA,1);      % spinodal (corrected)
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
    
    % find renormalized spinodal
    [chit(ii),phase(ii)]=spinodalRG(N,Nbar,FA);

    %%%%% LAM/HEX phase %%%%%
    if (phase(ii)==6 || phase(ii)==3)
%         chi13(ii)=lamhexRG(chis(ii),c(ii),d(ii),N,Nbar,miu,lam);
        chi13(ii)=chioot(chis(ii),c(ii),d(ii),N,Nbar,miu,lam,1,3);
    else
        chi13(ii)=0;
    end
    
    %%%%% HEX/BCC phase %%%%%
    if (phase(ii)==6)    
%         chi36(ii)=hexbccRG(chis(ii),c(ii),d(ii),N,Nbar,miu,lam);
        chi36(ii)=chioot(chis(ii),c(ii),d(ii),N,Nbar,miu,lam,3,6);
    else
        chi36(ii)=0;
    end
end
end

function chi13=chioot(chis,c,d,N,Nbar,miu,lam,n1,n2)
    % solver options
    options = optimset('Display','off',...
        'TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',1e14,'MaxIter',1e14);

    % initial guesses, lower bounds and upper bounds of solutions
    x02=[-1e2,1,1,1,1,1,1];
    lb2=[-1e5,0,0,0,0,0,0];
    ub2=[0,1e5,1e5,1e5,1e5,1e5,1e5];
    
    if ((n1==1 && n2==3) || (n1==3 && n2==1))
        % compare LAM/HEX phases
        theta1=0*miu;
        eta1=-lam/2;
        theta2=-miu;
        eta2=-lam/2;
    elseif ((n1==3 && n2==6) || (n1==6 && n2==3))
        % compare HEX/BCC phases
        theta1=-miu;
        eta1=-lam/2;
        theta2=-2*miu;
        eta2=3*lam/2;
    else
        error('Order-order transition between phases not defined')
    end

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

    F = @(x) [x(2)-x(1)-d*lam*power(x(2)*Nbar,-1/2),...
              x(5)-x(1)-d*lam*power(x(5)*Nbar,-1/2),...
              x(3)-x(1)-d*lam*power(x(3)*Nbar,-1/2)-n1*lam*x(4)^2,...
              x(6)-x(1)-d*lam*power(x(6)*Nbar,-1/2)-n2*lam*x(7)^2,...
              eta1*x(4)^2-theta1*x(4)+x(3),...
              eta2*x(7)^2-theta2*x(7)+x(6),...
              ((1/2/lam)*(x(3)^2-x(2)^2)+d*power(Nbar,-1/2)*(sqrt(x(3))-sqrt(x(2)))-...
                2*n1/3*theta1*x(4)^3+(1/2)*n1*eta1*x(4)^4)-...
              ((1/2/lam)*(x(6)^2-x(5)^2)+d*power(Nbar,-1/2)*(sqrt(x(6))-sqrt(x(5)))-...
                2*n2/3*theta2*x(7)^3+(1/2)*n2*eta2*x(7)^4)];
    [x,~] = lsqnonlin(F,x02,lb2,ub2,options);
    tau=x(1);
    chi13=chis-c^2*tau/(2*N);
end