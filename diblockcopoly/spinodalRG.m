function chi1=spinodalRG(N,NM,Nbar,FA)
% calculate the renormalized spinodals

% find spinodal
[chis,ks,~]=spinodal(N,NM,FA);

% find renormalized spinodal
chi1=spinodalfh(N,NM,Nbar,FA,ks,chis);
end

function chi1=spinodalfh(N,NM,Nbar,FA,ks,chis)
%% RG with no q dependence (F-H)

% calculate vertex functions
NQ=1;  % assume to Q dependence
[gam3,gam4]=calcgamma(N,NM,FA,NQ);
gam3=real(gam3);
gam4=real(gam4);

% calculate constant (estimate local second-order derivative)
xs=ks^2*(1/6)*N;
dx=xs*1e-3;
dks=dx/(2*ks*N/6);
d2c=N*(gamma2(N,NM,FA,ks+dks,0)-...
     2*gamma2(N,NM,FA,ks,0)+...
       gamma2(N,NM,FA,ks-dks,0))/(dx^2);
c=power(1/3*xs*d2c,1/2);

% parameters
miu = N*gam3/power(c,3);
lam = N*gam4(1)/power(c,4);
d = (3*xs)/(2*pi);

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

%% start solving self-consistent equations
fprintf('Step 3: Calculating renormalized spinodal at FA=%.2f, NM=%.2f\n',FA,NM)

% solver option
options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1e8,'MaxIter',1e8);

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
chi1=chis-c^2*tau/(2*N);
end