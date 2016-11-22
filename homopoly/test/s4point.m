%this code tests the calculation of 4-point correlation functions of wormlike chains
clear;
addpath(genpath('../chainstats'))

NMPLOT=1e1;     % Number of Kuhn steps
NMV=1e1;    %Number of Kuhn steps
if NMV>=1e2
%     QV = logspace(0,log10(50),51);
    QV = logspace(-2,2,501);
else
%     QV = logspace(0,log10(50),51);
    QV = logspace(-2,2,501);
end
% QV = QV/NMV(1);
% QV = QV/sqrt(r2wlc(NMV(1)));
TVtemp=linspace(0,pi/2,126);    %relative angles of wavevectors

% Plot options:
PLOT1 = 1;      % S vs q
PLOT2 = 1;      % Surface plot of S vs (q and angle)
readfolder = '../data/sdata';    % save to folder destination

% Calculation options:
CALCON = 0;     % Calculation on 
SAVEON = 0;     % Save to data
savefolder = '../data/sdata';    % save to folder destination

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
        s4temp=zeros(length(QV),length(TVtemp));

        for tt=1:length(TVtemp)
            T=TVtemp(tt);
            for ii=1:length(QV)
                Q=QV(ii);

                disp(['Calculating at N=',num2str(NM),', T=',num2str(T),' and Q=',num2str(Q)])

                %make wavevectors
                Q1=Q*[1,0,0];
                Q2=transpose(rotz(T)*Q1');
                Q3=-Q1;
                Q4=-Q2;

                %calculate s4
                s4temp(ii,tt)=s4wlc(NM,Q1,Q2,Q3,Q4,ORDEig,ORDL,ResLayer);
                s4temp(ii,tt)=real(s4temp(ii,tt)/power(NM,4));
            end
        end
        
        % post-process data according to symmetry (NOTE: only valid when
        % four vectors on same plane
        TV = linspace(0,2*pi,4*(length(TVtemp)-1)+1);
        s4=zeros(length(QV),length(TV));
        for tt=1:length(TV)
            T=TV(tt);
            for ii=1:length(QV)
                Q=QV(ii);
                
                IT = find(abs(TVtemp-T)<1e-2);
                if isempty(IT)
                    IT = find(abs(TVtemp-(pi-T))<1e-2);
                end
                
                if isempty(IT)
                    IT = find(abs(TVtemp-(T-pi))<1e-2);
                end
                
                if isempty(IT)
                    IT = find(abs(TVtemp-(2*pi-T))<1e-2);
                end
                
                if isempty(IT)
                    error('cannot apply symmetry')
                end
                IQ = find(QV==Q);
                s4(ii,tt) = s4temp(IQ,IT);
                
                if QV(ii)==0
                    s4(ii,tt) = 1;
                end
            end
        end
        
        if (SAVEON)
            foldername=strcat(savefolder,sprintf('/N%.2f',NM));
            if ~exist(foldername, 'dir')
                mkdir(foldername);
            end
            dlmwrite([foldername,'/s'],s4)
            dlmwrite([foldername,'/t'],TV)
            dlmwrite([foldername,'/k'],QV)
%             dlmwrite([foldername,'/k'],QV*NM)
%             dlmwrite([foldername,'/k'],QV*sqrt(r2wlc(NM)));
        end
    end
end
    
% %%%%%%%%%%%%%%% make plots %%%%%%%%%%%%%%%
% read in data
foldername=strcat(readfolder,sprintf('/N%.2f',NMPLOT));
s4=load([foldername,'/s']);
TV=load([foldername,'/t']);
QV=load([foldername,'/k']);

% optional: normalize s with kL
% for tt=1:length(TV)
%     s4(:,tt)=s4(:,tt).*QM'*NM;
% end

% plot 1: S vs q
if (PLOT1)
%     figure;plot(QV,s4(:,251)'.*QV,'k-')
    figure;plot(QV,s4(:,251)'.*QV*NMPLOT,'k-')
    set(gca,'xscale','log')
end

% plot 2: surface plot of S vs (q and angle)
if (PLOT2)
    figure;set(gca,'fontsize',50)
    surf(TV,QV,real(s4),'edgecolor','none','LineStyle','none','FaceLighting','phong');
%     ylabel('qL');xlabel('\theta');zlabel('S_{1234}')
%     ylabel('qR');xlabel('\theta');zlabel('S_{1234}')
    ylabel('Wavevector 2l_Pq');xlabel('Wavevector Angle \theta');zlabel('S_{1234}')
    set(gca,'xscale','linear');set(gca,'yscale','log');
    xlim([min(TV),max(TV)]);ylim([min(QV),max(QV)]);
    colorbar;view([0,90]);colormap(jet)
    set(gca, 'CLim', [0,1]);
    set(gca, 'XTick', [0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','2\pi/3','2\pi'})
    set(gca,'fontsize',20)
end