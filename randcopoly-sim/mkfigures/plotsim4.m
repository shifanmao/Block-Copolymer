function MU=plotsim4(EPS,LAM)
% Plots single chain correlations in Monte-Carlo simulation
% INPUTS::
%   EPS, number of Kuhn steps per monomer
%   LAM, degree of chemical correlation

% simulation folder
folder = '../../results/randcopoly-results/scalcbatch-12-15-15';
addpath('../functions/')

% simulation constants
M=8;        % number of blocks
G=5;  % number of discrete monomers

% range of chi params.
% chiind=fliplr([1:2:11,21,31,41]);
chiind = 1:41;
plotind=chiind(1:end);

% load simulation parameters
simparam=load([folder,'/chivals']);
% find corresponding simulation at EPS and LAM
SIMNUM = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
% finding simulation index(indices)
chemparam=load([folder,'/chemind']);
CHEMNUM=chemparam(chemparam(:,1)==SIMNUM,2);
% Flory-Huggins parameter
chiv = load(sprintf([folder,'/sdata-%d-%d/Sdata/chilist'],SIMNUM,CHEMNUM));
chiv = chiv/G;

folder = '../../results/randcopoly-results/MUdata';
cnt=find(chiind==plotind(1));
MU=zeros(length(plotind),1);
for ii = plotind
    CHI = chiv(ii);
    col = cnt/length(chiind);

    % Plot simulation results
    filename = sprintf([folder,'/MUMC_SIM%dCHEM%dCHI%.8f'],SIMNUM,CHEMNUM,CHI*G);
    data = load(filename);
    MU(cnt)=data;
    cnt = cnt+1;
end

figure;plot(chiv(plotind),MU,'o')
% savename = sprintf('../../results/randcopoly-results/random-simulation/pdata-lam%.2f.eps',LAM);
% saveas(gcf,savename,'epsc')