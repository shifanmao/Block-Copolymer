% plot radius of gyration of single chains
clear all;

% Plot Parameters
EPSV = [0.01,0.10,1.00];
LAM = 0.00;

col = ['r','b','k'];
f = figure(1);
hold;set(gca,'fontsize',30)
set(f,'position',[10,0,800,600])

for ii = 1:3
EPS = EPSV(ii);

% Simulation Parameters
N = 8;G = 5;Ree = 2;

% Calculate lengths of worm-like chains
L0 = Ree*EPS*((-0.5+0.5*exp(-2*EPS*G)+EPS*G)^(-0.5));  % Submonomer contour length

if EPS==0.01
    Lc = (N*G)*L0;  % Polymer contour length
else
    Lc = (N*G-1)*L0;  % Polymer contour length
end

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
% CHIV = linspace(0,4.5,length(R));

% Make a plot
if EPS == 5
    figure(1);plot(R(1:end-1,1),R(1:end-1,2)/sqrt(Rgid),'s-','linewidth',3,'color',col(ii))
    %figure(1);plot(CHIV,R(:,3)/sqrt(R2id),'s--','linewidth',3,'color',col(ii))
else
    figure(1);plot(R(:,1),R(:,2)/sqrt(Rgid),'s-','linewidth',3,'color',col(ii))
end
end

figure(1);xlabel('\chivG');ylabel('R_g/R_{id}')
ylim([.99,1.12])
set(gca,'YTickLabel',sprintf('%3.2f|',linspace(1,1.12,7)))
legend('N_M=0.05','N_M=0.50','N_M=5.00','location','northwest')
box on