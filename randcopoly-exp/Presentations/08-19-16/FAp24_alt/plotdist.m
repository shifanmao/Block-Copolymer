clear;close all

% OPTIONS
PLOTSIM = 1;    % plot simulation structure factor

% plot parameters
NREP = [1,80];  % number of replicas
NSNAP = 26:55;  % snapshots to average

% simulation parameters
EPS = 1.00;  % inter-bead segment rigidity (in unit of 2lp)

% simulation constants
boxl = 20;   % edge size of simulation
Ree = 2.0;   % average end-to-end distance of a monomer
G = 5;       % number of beads per monomer
N = 8;
FA = 0.5;    % chemical fraction of A species
NP = 2000;   % number of polymers
NK = (0:N*G-1)'/G;
L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);

% add/define paths
%loaddir = '../data/';    % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
savedir = 'savedata/';   % save directory
addpath('misc/');
addpath('../utility/');

% load CHI parameters
CHI = load(strcat(loaddir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

figure('Position', [100, 100, 1200, 900]);hold;set(gca,'fontsize',20);
figure('Position', [100, 100, 1200, 900]);hold;set(gca,'fontsize',20);
leg = {};

if (PLOTSIM)
    % plot simulation results
    for REP=NREP
      if max(NREP) == 1
        col = 1;
      else
        col = (REP-1)/(max(NREP)-1);
      end
      Ravg = [];Pavg = [];
      for SNAP=NSNAP
          fprintf('REP = %d, SNAP = %d\n',REP,SNAP)
          r=dlmread(strcat(loaddir,sprintf('r%dv%d',SNAP,REP)));

          % calculate bead-to-bead distance
          R = rcalc(r,NP);
          Ravg = [Ravg,R];

          % calculate end-to-end distribution
          [X,P]=pcalc(r,NP,L0);
          Pavg = [Pavg,P];
      end
      Ravg = mean(Ravg,2);
      Pavg = mean(Pavg,2);

      figure(1);plot(NK,sqrt(Ravg)/Ree,'o','linewidth',1.5,'color',[col 0 1-col]);
      figure(2);plot(X,Pavg,'o','linewidth',1.5,'color',[col 0 1-col]);
      leg{REP} = strcat('\chiG = ', sprintf('%.2f', CHI(REP)*G));
    end
end

figure(1);
% plot analytical solution
NK = logspace(-1,1,10);
NM = EPS*G;
DM = EPS*NK*G;
RM = NM-(1/2)*(1-exp(-2*NM));
RD = DM-(1/2)*(1-exp(-2*DM));
plot(NK,sqrt(RD/RM),'k-','linewidth',1.5)
plot([1e-1,1],[1,1],'k--','linewidth',1.5)
plot([1,1],[1e-1,1],'k--','linewidth',1.5)

% plot processing
set(gca,'xscale','log');set(gca,'yscale','log');box on
xlabel('J/G');ylabel('R_J/R_M')
legend(leg{NREP},'location','north');

figure(2);
pplot(N*G,EPS*N*G);
xlabel('R/L');ylabel('P(R/L)');box on
legend(leg{NREP},'location','north');
