%%%%%%%%%%%%%% DIAGRAM #1 %%%%%%%%%%%%%%
for EPS=[0.01,0.10,1.00]
% for EPS=[0.01]

phase = dlmread('simdata/phasediag','',1,0);
G = 5;
IND = find(phase(:,1)==EPS);
figure;hold;set(gca,'fontsize',20)

% LOCATIONS OF SNAPSHOTS
if EPS==0.01
    CHI1=[11.966,12.091,11.463,8.0645,5.1855];
    CHI2=[2.4505,4.0716,8.0645,10.805];
    LAM1=[-0.75:0.25:0.25];
    LAM2=ones(1,length(CHI2))*0;
    plot(LAM1,CHI1,'ko','linewidth',1.1,...
         'markersize',8,'markerfacecolor','k')
    plot(LAM2,CHI2,'k-s','linewidth',1.1,...
         'markersize',4)
elseif EPS==0.10
    CHI1=[13.049,12.889,11.666,8.0592,5.2472];
    LAM1=[-0.75:0.25:0.25];
    plot(LAM1,CHI1,'ko','linewidth',1.1,...
         'markersize',8,'markerfacecolor','k')
elseif EPS==1.00
    CHI1=[19.74,17.12,12.615,8.0044,5.1661];
    LAM1=[-0.75:0.25:0.25];
    plot(LAM1,CHI1,'ko','linewidth',1.1,...
         'markersize',8,'markerfacecolor','k')
end


% 
% % LOCATIONS OF SNAPSHOTS
% if EPS==0.01
%     CHI1=[11.236,11.966,12.091,11.463,8.0645,5.1855];
%     CHI2=[2.4505,4.0716,8.0645,10.805];
%     LAM1=[-1.00:0.25:0.25];
%     LAM2=ones(1,length(CHI2))*0;
%     plot(LAM1,CHI1,'ko','linewidth',1.1,...
%          'markersize',8,'markerfacecolor','k')
%     plot(LAM2,CHI2,'k-s','linewidth',1.1,...
%          'markersize',4)
% elseif EPS==0.10
%     CHI1=[13.08,13.049,12.889,11.666,8.0592,5.2472];
%     LAM1=[-1.00:0.25:0.25];
%     plot(LAM1,CHI1,'ko','linewidth',1.1,...
%          'markersize',8,'markerfacecolor','k')
% elseif EPS==1.00
%     CHI1=[22.068,19.74,17.12,12.615,8.0044,5.1661];
%     LAM1=[-1.00:0.25:0.25];
%     plot(LAM1,CHI1,'ko','linewidth',1.1,...
%          'markersize',8,'markerfacecolor','k')
% end

%LAM1=ones(1,length(CHI1))*-0.75;

% PLOT SECOND-ORDER ODT
if EPS == 0.01
    plot(phase(IND,2),phase(IND,4),'b-','linewidth',2);
    plot(phase(IND+3:IND(end),2),phase(IND+3:IND(end),4),'bo',...
        'markerfacecolor','b','linewidth',2);
else
    plot(phase(IND,2),phase(IND,4),'b-o',...
        'markerfacecolor','b','linewidth',2);
end
% PLOT FIRST-ORDER ODT
plot(phase(IND,2),phase(IND,5),'r-o',...
    'markerfacecolor','r','linewidth',2);
% PLOT MELTING CURVE
plot(phase(IND,2),phase(IND,8),'r--o',...
    'markerfacecolor','r','linewidth',2);
% PLOT PEAK EMERGENCE
plot(phase(IND,2),phase(IND,6),'m--o',...
    'markerfacecolor','r','linewidth',2);
% % PLOT MAX CHI
% plot(phase(IND,2),phase(IND,7),'k--',...
%     'markerfacecolor','r','linewidth',2);

% PLOT MF ODT (SECOND-ORDER)
filename = strcat('simdata/',...
    sprintf('MFspinNM%.2f',EPS*G));
MF = load(filename);
LAMV = MF(2:end,1);
CHIMF = MF(2:end,2);
plot(LAMV,CHIMF,'b:','linewidth',2);
% plot(LAMV,4*CHIMF,'b-','linewidth',2);

% INDICATE LIFSHITZ POINT
LAML = MF(1,1);
plot(LAML,0,'ks','markerfacecolor','k','markersize',8)
xpos = (LAML+1)/(1.5);
x = [xpos,xpos];y = ([0.18,0.15]);
% annotation('textarrow',x,y,'String',...
%     'Lifshitz Point','fontsize',18)

xlabel('Chemical Correlation \lambda')
ylabel('Flory-Huggins Parameter  \chiV_M')
% if EPS==0.01
    ylim([0,20]);
% else
%     ylim([0,25]);
% end
xlim([-1,.25]);box on
set(gca,'xtick',[-1.0,-0.75,-0.50,-0.25,0.00,0.25]);
set(gca,'xticklabels',{'-1.0','-0.75','-0.50','-0.25','0.00','0.25'})
set(gca,'linewidth',1.2)

savename = sprintf('PHASES_EPS%.2f.eps',EPS);
saveas(gcf,savename,'epsc')
% savename = sprintf('PHASES_EPS%.2f.png',EPS);
% saveas(gcf,savename,'png')
end
