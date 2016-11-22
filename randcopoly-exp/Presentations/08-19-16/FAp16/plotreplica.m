clear;
%close all

% plot options
NINTERV = 100;
NSTART = 1000;

% simulation parameters
G = 5;

% parse replicas
%loaddir = '../data/';  % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
addpath('misc/');
addpath('../utility/');

chi = load(strcat(loaddir,'chi'));
node = load(strcat(loaddir,'nodeNumber'));
[REPS,CHI,NTIME,NREP] = replica(chi,node);

figure;hold
for ii = 1:NREP
 plot(1:NINTERV:NTIME,CHI(1:NINTERV:NTIME,ii)*G,'linewidth',1.5)
end
xlabel('Monte Carlo Time');ylabel('\chiG');

figure;hold
for ii = NSTART:1:NTIME
    plot(CHI(NSTART,:)*G,CHI(ii,:)*G,'r.')
end
xlabel('\chi_0G');ylabel('\chi_hG')
