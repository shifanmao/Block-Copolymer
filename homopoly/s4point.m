%this code tests the calculation of 4-point correlation functions of wormlike chains
clear;
addpath('misc')

%Number of Kuhn steps
NMV=logspace(0,3,7);

%relative angles of wavevectors
TV=linspace(0,2*pi,1);

%Calculation parameters
ORDEig=4;
ORDL=6;
ResLayer=500;
ImagThreshold=1e-8;

for nn=1:length(NMV)
    NM=NMV(nn);
    
    %wavevector and structure factor
    QM=linspace(1e-5,10,50)/NM;
    s4=zeros(length(QM),length(TV));

    for tt=1:length(TV)
        T=TV(tt);
        for ii=1:length(QM)
            Q=QM(ii);
            
            disp(['Calculating at N=',num2str(NM),', T=',num2str(T),' and Q=',num2str(Q)])

            %make wavevectors
            Q1=Q*[1,0,0];
            Q2=transpose(rotz(T)*Q1');
            Q3=-Q1;
            Q4=-Q2;

            %calculate s4
            s4(ii,tt)=s4wlc(NM,Q1,Q2,Q3,Q4,ORDEig,ORDL,ResLayer);
            s4(ii,tt)=real(s4(ii,tt)/power(NM,4));
        end
    end

    foldername=sprintf('data/N%.2f',NM);
    if ~exist(foldername, 'dir')
        mkdir(foldername);
    end
    dlmwrite([foldername,'/s'],s4)
    dlmwrite([foldername,'/t'],TV)
    dlmwrite([foldername,'/k'],QM)
end

% % %%%%%%%%%%%%%%% make a surface plot %%%%%%%%%%%%%%%
% %Number of Kuhn steps
% NM=3.16;
% 
% % read in data
% foldername=sprintf('data/N%.2f',NM);
% s4=csvread([foldername,'/s']);
% TV=csvread([foldername,'/t']);
% QM=csvread([foldername,'/k']);
% 
% % optional: normalize s with kL
% for tt=1:length(TV)
%     s4(:,tt)=s4(:,tt).*QM';
% end
% 
% % plot
% figure;set(gca,'fontsize',15)
% surf(TV,QM,real(s4),'edgecolor','none','LineStyle','none','FaceLighting','phong');
% ylabel('kL');xlabel('\theta');zlabel('S_{1234}')
% set(gca,'xscale','linear');set(gca,'yscale','linear');
% xlim([min(TV),max(TV)]);ylim([min(QM),max(QM)]);
% set(gca, 'CLim', [0,2.]);colorbar;view([0,90])