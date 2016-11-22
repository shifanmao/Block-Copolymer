clear;
close all

global X1 Y1 X2 Y2 fitmin fitmax

% plot parameters
NEXP = [1:2:9];   % experiment plot range
%NREP = [19:1:24];  % number of replicas
NREP = [39:1:44];  % number of replicas

% simulation parameters
EPS = 1.00;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0.;     % degree of chemical correlation

% simulation constants
G = 5;       % number of beads per monomer
FA = 0.24;   % chemical fraction of A species

% experiment constants
TV = [22,40:20:180];  % temperature in degree C

% add/define paths
% %loaddir = '../data/';    % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
savedir = 'savedata/';   % save data directory
figdir = 'figures/';     % save figure directory

% define experiment data
sexp_full = load('expdata/PEG40MW1500.csv');

%% Figure 1: structure factors
CHI = load(strcat(loaddir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

% fitting wavevector q range
fitmin = 3e-2;
fitmax = 3e-1;

figure('Position', [100, 100, 1800, 900]);
hold;set(gca,'fontsize',20);leg = {};

X1 = [];X2 = [];
Y1 = [];Y2 = [];
cnt = 1;
for ii = NEXP
  if max(NEXP) == 1
    col = 1;
  else
    col = (ii-1)/(max(NEXP)-1);
  end

  % load experiment data
  x1 = sexp_full(:,1);
  y1 = sexp_full(:,ii+1);

  % load simulation results
  REP = NREP(cnt);
  SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
  ssim = load(strcat(savedir,SAVEFILENAME));

  % log distributed number of points
  INDSIM = unique(round(logspace(0,log10(length(ssim)),50)));
  x2 = ssim(INDSIM,1);
  y2 = ssim(INDSIM,2);

  % append to matrix
  X1 = [X1,x1];Y1 = [Y1,y1];
  X2 = [X2,x2];Y2 = [Y2,y2];

  cnt = cnt+1;
end

% use non-linear fitting
RM0 = 20; SCALE0 = 0.5; BG0 = ones(1,length(NEXP))*3.3;
x0 = [RM0,SCALE0,BG0];
lb = [1,1e-2,ones(1,length(NEXP))*0];
ub = [50,100,ones(1,length(NEXP))*100];
options = optimset('MaxFunEvals',500,'Display','off','TolFun',1e-8);
x = lsqnonlin(@fits_rm,x0,lb,ub,options);

% extract fitting results
RM = x(1)
SCALE = x(2)
BG = x(3:end)

cnt = 1;
for ii = NEXP
  if max(NEXP) == 1
    col = 1;
  else
    col = (ii-1)/(max(NEXP)-1);
  end

  % plot experiment data
  x1 = sexp_full(:,1);
  y1 = sexp_full(:,ii+1);
  plot(x1*RM,y1,'linewidth',2,'color',[col 0 1-col])

  leg{cnt} = strcat('T = ', sprintf('%dC^o', TV(ii)));
  cnt = cnt+1;
end

cnt = 1;
for ii = NEXP
  if max(NEXP) == 1
    col = 1;
  else
    col = (ii-1)/(max(NEXP)-1);
  end

  REP = NREP(cnt);
  SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
  ssim = load(strcat(savedir,SAVEFILENAME));
  x2 = ssim(1:end-1,1);
  y2 = ssim(1:end-1,2);

  plot(x2,exp(SCALE*log(y2)+BG(cnt)),'.-','MarkerSize',15,...
       'linewidth',1.5,'color',[col 0 1-col])

  leg{cnt+length(NEXP)} = strcat('\chiG = ', sprintf('%.2f', CHI(REP)*G));
  cnt = cnt+1;
end

cnt = 1;
for ii = NEXP
  if max(NEXP) == 1
    col = 1;
  else
    col = (ii-1)/(max(NEXP)-1);
  end

  % plot experiment data
  x1 = sexp_full(:,1);
  y1 = sexp_full(:,ii+1);
  
  % find out fitting range
  xmin = fitmin; xmax = fitmax;
  minind = find(x1>=xmin);minind = minind(1);
  maxind = find(x1<=xmax);maxind = maxind(end);

  plot(x1(minind)*RM,y1(minind),'k.','MarkerSize',30);
  plot(x1(maxind)*RM,y1(maxind),'k.','MarkerSize',30);

  cnt = cnt+1;
end


% plot process
%legend(leg{NEXP});
legend(leg);
axis([0.3,10,5e0,5e2])
txt = 'fit q_0';
%text(0.65,16,txt,'fontsize',15);
xlabel('R_Mq');ylabel('S(a.u.)');box on
set(gca,'xscale','log');set(gca,'yscale','log')
saveas(gcf,strcat(figdir,'sk'),'epsc')
