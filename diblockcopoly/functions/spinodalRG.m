function [chit,phase]=spinodalRG(N,Nbar,FA)
% Find the renormalized spinodal by FH theory
% Usage: chi1=spinodalRG(N,Nbar,FA)

% find spinodal
[chis,ks,~]=spinodal(N,FA);

% find renormalized spinodal of each phase
chi1=spinodalfh(N,Nbar,FA,ks,chis,1);
chi3=spinodalfh(N,Nbar,FA,ks,chis,3);
chi6=spinodalfh(N,Nbar,FA,ks,chis,6);

% find renormalized spinodal
chiall=[chi1,chi3,chi6];
phases=[1,3,6];
[chi,order]=sort(chiall);
chit=chi(1);
phase=phases(order(1));
end

function chit=spinodalfh(N,Nbar,FA,ks,chis,n)
% calculate renormalized spinodal from FH theory
    % calculate vertex functions
    NQ=1;  % assume to Q dependence
    [gam3,gam4]=calcgamma(N,FA,NQ);
    gam3=real(gam3);
    gam4=real(gam4);
    
    % calculate constant (estimate local second-order derivative)
    xs=ks^2*(1/6)*N;
    dx=xs*1e-3;
    dks=dx/(2*ks*N/6);
    d2c=N*(gamma2(N,FA,ks+dks,0)-2*gamma2(N,FA,ks,0)+gamma2(N,FA,ks-dks,0))/(dx^2);
	c=power(1/3*xs*d2c,1/2);
    
    % parameters
    miu = N*gam3/power(c,3);
    lam = N*gam4(1)/power(c,4);
    d = (3*xs)/(2*pi);
    
    %%%%% LAM/HEX/BCC phase %%%%%
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

    %% start solving self-consistent equations
    fprintf('Step 2: Calculating renormalized spinodal at FA=%.2f, N=%.2e\n',FA,N)
    
    % solver option
    options = optimset('Display','off',...
        'TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',1e12,'MaxIter',1e12);

    % initial guesses, lower bounds and upper bounds of solutions
    x0 = [-1e2,1,1,1];
    lb=[-1e7,0,0,0];
    ub=[0,1e5,1e5,1e5];
    
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

    F = @(x) [x(2)-x(1)-d*lam*power(x(2)*Nbar,-1/2),...
              x(3)-x(1)-d*lam*power(x(3)*Nbar,-1/2)-n*lam*x(4)^2,...
              eta*x(4)^2-theta*x(4)+x(3),...
              (1/2/lam)*(x(3)^2-x(2)^2)+...
              d*power(Nbar,-1/2)*(sqrt(x(3))...
              -sqrt(x(2)))-2*n/3*theta*x(4)^3+(1/2)*n*eta*x(4)^4];
    [x,~] = lsqnonlin(F,x0,lb,ub,options);
    tau=x(1);
    chit=chis-c^2*tau/(2*N);
end