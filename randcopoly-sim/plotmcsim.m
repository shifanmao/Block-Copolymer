function [chis,chiplot,ks,ksim,sinv_theory,sinv_sim]=plotmcsim(EPS,LAM,PLOTON)
%% Plots density-density correlation of sim. and theory
% INPUTS::
%   EPS = number of Kuhn steps per monomer
%   LAM = degree of chemical correlation

% simulation folder
% folder = 'scalcbatch-09-30-15';
folder = 'scalcbatch-12-15-15';

% simulation constants
FA=0.5;  % fraction of A blocks
M=8;        % number of blocks
G=5;  % number of discrete monomers
NM=G*EPS;  % number of Kuhn steps per monomer

% range of chi params.
chiind=fliplr([1:2:11,22,32,42]);
plotind=chiind(1:end);
% chiind=fliplr(1:46);
% plotind=chiind(5:end);

% load simulation parameters
simparam=load([folder,'/chivals']);
% find corresponding simulation at EPS and LAM
ind = find(simparam(:,1)==EPS & simparam(:,2)==LAM);
% finding simulation index(indices)
chemparam=load([folder,'/chemind']);
SIMNUM = ind;
CHEMNUM=chemparam(chemparam(:,1)==ind,2);
% Flory-Huggins parameter
CHIV = load(sprintf([folder,'/sdata-%d-%d/Sdata/chilist'],SIMNUM,CHEMNUM));
CHIV = CHIV/G;

% Find spinodal CHI and critical k
[chis,ks,~]=spinodal(M,NM,LAM,FA);
R2=-0.5+0.5*exp(-2*NM)+NM;
k=logspace(log10(0.6),log10(12),100)./sqrt(R2); % wavevectors

% %% Figure 1: density-density correlation
if PLOTON==1
    f1=figure;hold;set(gca,'fontsize',30)
    set(f1,'position',[0,0,800,600])
    cnt=find(chiind==plotind(1));
    for ii = plotind
        CHI = CHIV(ii);
        col = cnt/length(chiind);

        if CHI<chis*EPS*0.85
            % Plot analytical theory (note conversion of monomer volume v, factor of EPS)
            val=gamq2wlc(M,NM,LAM,k);
            plot(k.*sqrt(R2),1./(-2*CHI+EPS*NM*val/(FA*(1-FA))),'--','color',[1-col 0 col],'linewidth',3)

            % Plot power law scaling
            % if cnt==length(chiind);
            %     plotpower(-1,k.*sqrt(r2),1./(G*params(ind,1)-2*CHI/g),5,8,0.85)
            % end
        end

        % Plot simulation results
        filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
            SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
        S = load(filename);
        plot(S(:,1),S(:,2),'.-','color',[1-col 0 col],'linewidth',3,'markersize',20)
        cnt = cnt+1;
    end

    xlabel('R_Mq');
    ylabel('<\psi(q)\psi(-q)>')
    % xlabel('Wavevector \it{qR_M}');
    % ylabel('Density correlation \it{<\psi(q)\psi(-q)>}')
    xlim([.6,12]);ylim([1e-1,2e2]);
    set(gca,'Ytick',[1e-1 1e0 1e1 1e2])
    set(gca,'Xtick',[1,10])
    % plot([ks,ks]*sqrt(R2),[1e-1,2e2],'k--');
    
    set(gca,'xscale','log');set(gca,'yscale','log')
    box on
end

%% Figure 2: S^{-1}(q*)
% range of chi params.
chiind=1:46;
plotind=1:2:46;

ksim = zeros(length(plotind),1);
sinv_theory = zeros(length(plotind),1);
sinv_sim = zeros(length(plotind),1);
chiplot = CHIV(plotind)*G;
cnt = 1;
for ii = plotind
    CHI = CHIV(ii);
    
%    if CHI<=chis*EPS
        % Plot analytical theory (note conversion of monomer volume v, factor of EPS)
        val=gamq2wlc(M,NM,LAM,ks);
        sinv_theory(cnt)=-2*CHI+EPS*NM*val/(FA*(1-FA));
%    end

    % Plot simulation results
    filename = sprintf([folder,'/sdata-%d-%d/Sdata/SMC_SIM%dCHEM%dCHI%.8f'],...
        SIMNUM,CHEMNUM,SIMNUM,CHEMNUM,CHI*G);
    S = load(filename);
    
    % Find peak position
    ind = find(S(:,2)==max(S(:,2)));
    ksim(cnt) = S(ind(1),1);
    sinv_sim(cnt) = 1./S(ind(1),2);
    cnt = cnt+1;
    
    savename = sprintf('structure-figures/sfig-eps%.2f-lam%.2f.eps',EPS,LAM);
    saveas(gcf,savename,'epsc')
end