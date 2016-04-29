function plotsim3(EPS,LAM,PLOTON)
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
L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);
LP=L0/(2*EPS);
N=M*G;
RM=2;
NK=[0:N-1]'/G;

% range of chi params.
% chiind=fliplr([1:2:11,21,31,41]);
% plotind=chiind(1:end);
chiind = [41,1];
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

folder = '../../results/randcopoly-results/Rdata';
cnt=find(chiind==plotind(1));

for ii = plotind
    CHI = chiv(ii);
    
    EPSV = [0.01,0.10,1.00];
    coli = find(EPSV==EPS);
    col = (coli-1)/(length(EPSV)-1);

    % Plot simulation results
    filename = sprintf([folder,'/RMC_SIM%dCHEM%dCHI%.8f'],SIMNUM,CHEMNUM,CHI*G);
    RB = load(filename);
    
    if PLOTON
        if ii==1
            %plot(NK,sqrt(RB)/RM,'-+','markersize',6,'color',[col 0 1-col],'linewidth',2)
            plot(NK,sqrt(RB)/RM,'+','markersize',6,'color',[col 0 1-col])
        else
            %plot(NK,sqrt(RB)/RM,'-o','markersize',6,'color',[col 0 1-col],'linewidth',2)
            plot(NK,sqrt(RB)/RM,'o','markersize',6,'color',[col 0 1-col])
        end
    end
    
    cnt = cnt+1;
end

% % Plot theoretical results
% R2 = NK-(1/2)*(1-exp(-2*NK));

% savename = sprintf('../../results/randcopoly-results/random-simulation/pdata-lam%.2f.eps',LAM);
% saveas(gcf,savename,'epsc')