%this code tests the calculation of 3-point vertex functions of wormlike chains
clear;
addpath(genpath('../chainstats'))

% Plot options:
PLOT1 = 1;      % S vs q
PLOT2 = 1;      % Surface plot of S vs (q and angle)
NMPLOT=1e4;     % Number of Kuhn steps
readfolder = '../data/gamdata';

% Calculation options:
CALCON = 0;     % Calculation on 
SAVEON = 0;     % Save to data
NM_WLC=1e4;     % Number of Kuhn steps
savefolder = '../data/gamdata';
QM=logspace(-2,2,50)';
ang = linspace(0,2*pi,51);

if (CALCON)
    gam3 = gamma3(QM,ang,NM_WLC);
    
    if (SAVEON)
        if (length(ang) == 1)
            filename = strcat(savefolder,sprintf('/gam3N1e%d',log10(NM_WLC)));
            dlmwrite(filename,[log10(QM),log10(gam3)],'precision','%.2f')
        else
            filename = strcat(savefolder,sprintf('/gam3N1e%dang',log10(NM_WLC)));
            data = horzcat([0;log10(QM)],vertcat(ang,real(log10(gam3))));
            dlmwrite(filename,data,'precision','%.2f')
        end
    end
end

if (PLOT1)
    % load option
    filename = strcat(readfolder,sprintf('/gam3N1e%d',log10(NMPLOT)));
    x = load(filename);
    QM = 10.^(x(:,1));
    gam3 = 10.^(x(:,2));

    figure('Position', [100, 100, 1200, 900])
    hold;set(gca,'fontsize',50)
    plot(QM,gam3,'b-','linewidth',6);
    plot(QM,4./power(QM*sqrt(NM_WLC/6),4),'k--','linewidth',4);
    set(gca,'xscale','log');set(gca,'yscale','log');box on 
    xlabel('Wavevector K');ylabel('Cubic-order Vertex \Gamma^{(3)}')

    % save to figure
    filename = strcat(savefolder,sprintf('/gam3N1e%d',log10(NM_WLC)));
    saveas(gcf,strcat(filename,'.eps'),'epsc')
end

if (PLOT2)
    % load option
    filename = strcat(readfolder,sprintf('/gam3N1e%dang',log10(NMPLOT)));
    if exist(filename,'file')
        x = load(filename);
        QM = 10.^(x(2:end,1));
        ang = x(1,2:end);
        gam3 = 10.^(x(2:end,2:end));

        figure('Position', [100, 100, 1200, 900])
        hold;set(gca,'fontsize',50)
        surf(ang,QM,gam3);colormap(jet)
        set(gca,'yscale','log');set(gca,'zscale','log')
        xlim([0,2*pi]);view([0,90])
        xlabel('\theta');ylabel('Wavevector K');
        title('Cubic-order Vertex \Gamma^{(3)}')

        % save to figure
        filename = strcat(savefolder,sprintf('/gam3N1e%dang',log10(NM_WLC)));
        saveas(gcf,strcat(filename,'.eps'),'epsc')
    else
        warning('Cannot make surface plot, since file non-existent')
    end
end