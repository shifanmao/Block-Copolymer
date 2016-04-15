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
NK=[1:N-1]*EPS;

% range of chi params.
% chiind=fliplr([1:2:11,21,31,41]);
% plotind=chiind(1:end);
chiind = [1,41];
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

if PLOTON
%     fig1=figure;hold;set(gca,'fontsize',18)
%     fig2=figure;hold;set(gca,'fontsize',18)
%     figure(1);hold;
%     figure(2);hold;
end
for ii = plotind
    CHI = chiv(ii);
    col = (cnt-1)/(length(chiind)-1);

    % Plot simulation results
    filename = sprintf([folder,'/RMC_SIM%dCHEM%dCHI%.8f'],SIMNUM,CHEMNUM,CHI*G);
    RB = load(filename);
    RB = RB(2:end);
    
    if PLOTON
        %figure(fig1);
        figure(1);
        %plot(1:N,RB/((L0).^2),'-o','color',[col 0 1-col],'linewidth',2,'markersize',5)
        EPSV=[0.01,0.10,1.00];
        icol = find(EPSV==EPS);
        col = (icol-1)/(length(EPSV)-1);
        if ii==1
            plot(1:N-1,sqrt(RB/((L0).^2)),'-x','markersize',5,'color',[col 0 1-col],'linewidth',1)
        else
            plot(1:N-1,sqrt(RB/((L0).^2)),'-o','markersize',5,'color',[col 0 1-col],'linewidth',1)
        end

        %figure(fig2);
        figure(2);
        if ii==1
            plot(NK,sqrt(RB/((2*LP).^2)),'-x','markersize',5,'color',[col 0 1-col],'linewidth',1)
        else
            plot(NK,sqrt(RB/((2*LP).^2)),'-o','markersize',5,'color',[col 0 1-col],'linewidth',1)
        end
    end
    
    cnt = cnt+1;
end

% Plot theoretical results
NK=[1:N]*EPS;
R2 = NK-(1/2)*(1-exp(-2*NK));
if PLOTON
    % figure(fig1);
    figure(1);
%     plot(1:N,power(1:N,2),'k-')

    % figure(fig2);
    figure(2);
%     plot(NK,R2,'k-','linewidth',2)


end

% savename = sprintf('../../results/randcopoly-results/random-simulation/pdata-lam%.2f.eps',LAM);
% saveas(gcf,savename,'epsc')