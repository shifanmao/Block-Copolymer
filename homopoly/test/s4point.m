%this code tests the calculation of 4-point correlation functions of wormlike chains
clear;
addpath(genpath('../chainstats'))

% Plot options:
PLOT1 = 1;      % S vs q
PLOT2 = 1;      % Surface plot of S vs (q and angle)
NMPLOT=1e4;     % Number of Kuhn steps
readfolder = '../data';    % save to folder destination

% Calculation options:
CALCON = 1;     % Calculation on 
SAVEON = 1;     % Save to data
savefolder = '../data';    % save to folder destination
NMV=1e4;    %Number of Kuhn steps
% NMV=logspace(0,3,7);
QV=linspace(0,500,101); %wavevectors in unit of L contour length
TV=linspace(0,2*pi,101);    %relative angles of wavevectors

% %%%%%%%%%%%%%%% start calculations %%%%%%%%%%%%%%%
if (CALCON)
    %Calculation parameters
    ORDEig=2;
    ORDL=2;
    ResLayer=500;
    ImagThreshold=1e-8;

    for nn=1:length(NMV)
        NM=NMV(nn);

        %wavevector and structure factor
        QM=QV/NM;
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

        if (SAVEON)
            foldername=strcat(savefolder,sprintf('/N%.2f',NM));
            if ~exist(foldername, 'dir')
                mkdir(foldername);
            end
            dlmwrite([foldername,'/s'],s4)
            dlmwrite([foldername,'/t'],TV)
            dlmwrite([foldername,'/k'],QM)
        end
    end
end
    
% %%%%%%%%%%%%%%% make plots %%%%%%%%%%%%%%%
% read in data
foldername=strcat(readfolder,sprintf('/N%.2f',NMPLOT));
s4=load([foldername,'/s']);
TV=load([foldername,'/t']);
QM=load([foldername,'/k']);

% optional: normalize s with kL
% for tt=1:length(TV)
%     s4(:,tt)=s4(:,tt).*QM'*NM;
% end

% plot 1: S vs q
if (PLOT1)
    figure;plot(QM*NMPLOT,s4,'k--')
end

% plot 2: surface plot of S vs (q and angle)
if (PLOT2)
    figure;set(gca,'fontsize',15)
    surf(TV,QM,real(s4),'edgecolor','none','LineStyle','none','FaceLighting','phong');
    ylabel('kL');xlabel('\theta');zlabel('S_{1234}')
    set(gca,'xscale','linear');set(gca,'yscale','linear');
    xlim([min(TV),max(TV)]);ylim([min(QM),max(QM)]);
    colorbar;view([0,90]);colormap(jet)
    % set(gca, 'CLim', [0,2.]);
end