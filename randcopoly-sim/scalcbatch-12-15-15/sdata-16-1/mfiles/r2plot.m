% plot radius of gyration of single chains
clear all;

% Plot Parameters
EPSV = [0.01,0.10,1.00];
LAM = -0.25;

col = ['r','b','k'];
figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',20)

for ii = 1:3
EPS = EPSV(ii);

% Simulation Parameters
N = 8;G = 5;Ree = 2;

% Calculate lengths of worm-like chains
L0 = Ree*EPS*((-0.5+0.5*exp(-2*EPS*G)+EPS*G)^(-0.5));  % Submonomer contour length
Lc = (N*G-1)*L0;  % Polymer contour length
Lp = L0/2/EPS;  % Persistent length
Rgid = (1/3)*Lp*Lc-Lp^2+2*Lp^3/Lc*(1-(Lp/Lc)*(1-exp(-Lc/Lp)));

Lm = (G-1)*L0;  % Monomer contour length
%R2id = (1/3)*Lp*Lm-Lp^2+2*Lp^3/Lm*(1-(Lp/Lm)*(1-exp(-Lm/Lp)));
%R2id = (2*Lp)^2*(Lm/2/Lp-0.5*(1-exp(-Lm/Lp)));
%R2id = sqrt(R2id)

% Load data
simparam=load('chivals');
% find corresponding simulation at EPS and LAM
ind = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
% finding simulation index(indices)
chemparam=load('chemind');
SIMNUM = ind;
CHEMNUM=chemparam(chemparam(:,1)==ind,2);

% Read in data
filename = sprintf('rdata/RMC_SIM%dCHEM%d',SIMNUM,CHEMNUM)
R = load(filename);
CHIV = linspace(0,4.5,length(R));

% Make a plot
figure(1);plot(CHIV,R(:,2)/sqrt(Rgid),'s-','linewidth',3,'color',col(ii))
%figure(1);plot(CHIV,R(:,3)/sqrt(R2id),'s--','linewidth',3,'color',col(ii))
end

figure(1);xlabel('\chi/\chi_S');ylabel('R_g/R_{id}')
