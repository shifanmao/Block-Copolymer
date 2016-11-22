clear;
%close all

% plot parameters
NEXP = [1,3];   % experiment plot range
NREP = [39,42];  % number of replicas

% simulation parameters
EPS = 1.00;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0.;     % degree of chemical correlation

% simulation constants
G = 5;       % number of beads per monomer
FA = 0.24;   % chemical fraction of A species

% experiment constants
TV = [22,40:20:180];  % temperature in degree C
rm=32.3348;

% add/define paths
%loaddir = '../data/';    % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
savedir = 'savedata/';   % save data directory
figdir = 'figures/';     % save figure directory

%% Figure 1: structure factors
CHI = load(strcat(loaddir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

% load experiment data
sexp_full = load('expdata/PEG40MW1500.csv');

figure('Position', [100, 100, 1200, 900]);
hold;set(gca,'fontsize',20);leg = {};

Y1 = [];YP = [];
cnt = 1;
for ii = NEXP
  if max(NEXP) == 1
    col = 1;
  else
    col = (ii-1)/(max(NEXP)-1);
  end
  x1 = sexp_full(:,1)*rm;
  y1 = sexp_full(:,ii+1);
  plot(x1,y1,'linewidth',2,'color',[col 0 1-col])

  REP = NREP(cnt);
  SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
  ssim = load(strcat(savedir,SAVEFILENAME));
  x2 = ssim(1:end-1,1);
  y2 = ssim(1:end-1,2);

  % start fitting adjustment parameters
  xmin = max([1.4,x2(1)]);minind = find(x1>=xmin);minind = minind(1);
  xmax = min([10,x2(end)]);maxind = find(x1<=xmax);maxind = maxind(end);
  x1 = x1(minind:maxind);
  y1 = y1(minind:maxind);
  x1star(cnt) = x1(find(y1 == max(y1)));  % find peak
  yp = [ones(length(x1),1),interp1(x2,log(y2),x1)];

  Y1 = [Y1;log(y1)];
  YP = [YP;yp];

  leg{cnt} = strcat('T = ', sprintf('%dC^o', TV(ii)));cnt = cnt+1;
end

% use linear regression
B = YP \ Y1;

% calculate residual
Ysim = B(1)+B(2)*YP(:,2);
Yresid = Y1 - Ysim;
SSresid = sum(Yresid.^2)

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
  x2star(cnt) = x2(find(y2 == max(y2)));  % find peak

%  plot(x2,B(1)+B(2)*y2,'.-','MarkerSize',15,'linewidth',1.5,'color',[col 0 1-col])
  plot(x2,exp(B(1)+B(2)*log(y2)),'.-','MarkerSize',15,'linewidth',1.5,'color',[col 0 1-col])
  leg{cnt+length(NEXP)} = strcat('\chiG = ', sprintf('%.2f', CHI(REP)*G));cnt = cnt+1;
end

% plot process
%legend(leg{NEXP});
legend(leg);
axis([0.5,15,8e0,5e2])
xlabel('R_Mq');ylabel('S(a.u.)');box on
set(gca,'xscale','log');set(gca,'yscale','log')
saveas(gcf,strcat(figdir,'sk'),'epsc')
