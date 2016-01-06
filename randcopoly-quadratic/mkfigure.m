clear;
close all

% input files
WLCfile = 'data/WLC_N100NNM11.mat';
GSfile = 'data/GS_N100NM1e6';
RRfile = 'data/RR_N100NM1.mat';

% Range of NM
load(WLCfile)
ind = [1:length(NMV)];

%% plot 1a (spinodal vs. LAM)
f2a=figure;hold;set(gca,'fontsize',20)
set(f2a,'position',[0,0,800,600])
% plot worm like chain
load(WLCfile)
NMR=NMV(ind);
for NM = NMR
    ii = find(NMV==NM);
    col = (ii-1)/(length(ind)-1);
    plot(LAMV,CHIS(:,ii),'color',[col 0 1-col],'linewidth',3.)
end

% plot Rigid rod
load(RRfile);
plot(LAMV,CHIS,'k--','linewidth',3.);

% plot Gaussian chain
load(GSfile);
plot(LAMV,CHIS,'k-','linewidth',3.);

xlabel('Chemical Correlation \it{\lambda}');
ylabel('Spinodal 4f_Af_Bv\chi_SN_m')
box on

%% plot 2b (ks vs LAM)
f2b = figure;hold;set(gca,'fontsize',20)
set(f2b,'position',[100,0,800,600])
% plot worm like chain
load(WLCfile)
NMR=NMV(ind);
for NM = NMR
    ii = find(NMV==NM);
    col = (ii-1)/(length(ind)-1);
    plot(LAMV,sqrt(R2(ii)).*ks(:,ii),'color',[col 0 1-col],'linewidth',3.)
end

% plot Rigid rod
load(RRfile);
plot(LAMV,100*ks,'k--','linewidth',3.);

% plot Gaussian chain
load(GSfile);
plot(LAMV,sqrt(1e6)*ks,'k-','linewidth',3.);

xlabel('Chemical Correlation \lambda');
ylabel('Critical wavevector R_mq^*')
box on
h2 = gca;
yticks=0:0.5:4.5;
set(h2(1),'YTick',yticks)
set(h2(1),'YTickLabel',sprintf('%.1f|',yticks))

%% plot 2b inset (D/R_M vs LAM)
pos = 0.7;
h = axes('position',ones(4,1)*pos);hold
set(gca,'fontsize',18)
% plot worm like chain
load(WLCfile)
NMR=NMV(ind);
for NM = NMR
    ii = find(NMV==NM);
    col = (ii-1)/(length(ind)-1);
    plot(LAMV,2*pi./(sqrt(R2(ii)).*ks(:,ii)),'color',[col 0 1-col],'linewidth',3.)
end

% plot Rigid rod
load(RRfile);
plot(LAMV,2*pi./(100*ks),'k--','linewidth',3.);

% plot Gaussian chain
load(GSfile);
plot(LAMV,2*pi./(sqrt(1e6)*ks),'k-','linewidth',3.);

p1=0.55;p2=0.55;pend=0.88;
set(h,'position',[p1,p2,pend-p1,pend-p2])
% xlabel('Chemical Correlation \lambda');
% ylabel('Critical wavelength D/R_M')
xlabel('\it{\lambda}');
ylabel('\it{D^*/R_M}')
ylim([1.2,3]);xlim([-1,0]);box on;
yticks=[1.5,2.0,2.5,3.0];
set(h(1),'YTick',yticks)
set(h(1),'YTickLabel',sprintf('%.1f|',yticks))

%% plot 6 susceptibility vs LAM/NM
f2a=figure;hold;set(gca,'fontsize',20)
set(f2a,'position',[0,0,800,600])
% plot worm like chain
load(WLCfile)
NMR=NMV(ind);
for NM = NMR
    ii = find(NMV==NM);
    col = (ii-1)/(length(ind)-1);
    plot(LAMV,log(-d2gam2(:,ii)),'color',[col 0 1-col],'linewidth',3.)
end

% plot Rigid Rod
load(RRfile);
plot(LAMV,log(-d2gam2),'k--','linewidth',3.);

% plot Gaussian chain
load(GSfile);
plot(LAMV,log(-d2gam2),'k-','linewidth',3.);

xlabel('Chemical Correlation \it{\lambda}');
ylabel('log(|Susceptibility|)')
box on