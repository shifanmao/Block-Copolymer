clear;close all

% OPTIONS
PLOTSIM = 1;    % plot simulation structure factor
PLOTMF = 0;     % plot mean-field structure factor
SAVESIM = 1;    % save simulation structure factor to file

% plot parameters
NREP = [1:80];  % number of replicas
NSNAP = 26:100;  % snapshots to average
lksample = 10;

% simulation parameters
EPS = 1.00;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0.;     % degree of chemical correlation

% simulation constants
boxl = 20;   % edge size of simulation
Ree = 2.0;   % average end-to-end distance of a monomer
G = 5;       % number of beads per monomer
FA = 0.24;   % chemical fraction of A species

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

%% (Structure factors)
if (PLOTSIM || PLOTMF)
  figure;hold;set(gca,'fontsize',20);
end

if (PLOTSIM)
    % plot simulation results
    for REP=NREP
      if max(NREP) == 1
        col = 1;
      else
        col = (REP-1)/(max(NREP)-1);
      end
      savg = [];
      for SNAP=NSNAP
          fprintf('REP = %d, SNAP = %d\n',REP,SNAP)
          r=dlmread(strcat(loaddir,sprintf('r%dv%d',SNAP,REP)));
          [k,s]=scalc(r,boxl,lksample);
          savg = [savg,s];
      end
      K = k*Ree; S = mean(savg,2);
      plot(K,S,'o-','linewidth',1.5,'color',[col 0 1-col]);

      % save to file
      if (SAVESIM)
          SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
          dlmwrite(strcat(savedir,SAVEFILENAME),[K,S]);
      end
        
      % fit to Lorentzian
      %[KS,SINV,D2S,ERR]=calcserr(K,S,G);
        
      % find peak location
      %[pks,locs] = findpeaks(S);
    end
end

if (PLOTMF)
    % Plot mean-field theory
    % load MF results
    filename = sprintf('../utility/sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
    MF = load(filename);
    CHIS = MF(1,1);   %spinodal
    for REP=NREP
        if max(NREP)==1
            col = 1;
        else
            col = (REP-1)/(max(NREP)-1);
        end

        if CHI(REP) < CHIS/G
            K = MF(2:end,1); %wavevectors
            S = 1./(-2*CHI(REP)+EPS*MF(2:end,2));
            plot(K,S,'--','color',[col 0 1-col],'linewidth',2);
        end
    end
end

% plot processing
xlabel('R_Mq');ylabel('S(q)');box on
set(gca,'xscale','log');set(gca,'yscale','log')
