clear;
% close all

CALCON = 0;

% universal parameters
N = 100;

filename1 = '../data/figureS5A';
filename2 = '../data/figureS5B';
filename3 = '../data/figureS5C';
filename4 = '../data/figureS5D';
filename5 = '../data/figureS5E';
filename6 = '../data/figureS5F';
filenames = {filename1, filename2, filename3, filename4, filename5, filename6};

if CALCON
    FAV = linspace(.1, .9, 11);
    LAMV = linspace(-1, 1, 11);
    CHIABSV = zeros(length(FAV), length(LAMV));
    KSV = zeros(length(FAV), length(LAMV));

    cnt = 1;
    for PHIP = [.2, .5, .8]
        for NM = [1e-2, 1e2]
            for ii = 1:length(FAV)
                FA = FAV(ii)
                for jj = 1:length(LAMV)
                    LAM = LAMV(jj);

                    [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
                    CHIABSV(ii, jj) = CHIABS*NM*PHIP;
                    KSV(ii, jj) = KS*sqrt(r2(NM));
                end
            end

            filename = filenames{cnt};
            save(filename, 'CHIABSV', 'KSV', 'LAMV', 'FAV', 'PHIP', 'NM');
            cnt = cnt+1;
        end
    end
    
else
    for ff = 1:2
        load(filenames{ff})

        for ii = 1:length(FAV)
            FA = FAV(ii);
            LAMMIN = max([-(1-FA)/FA, -FA/(1-FA)]);
            IND = find(LAMV<LAMMIN);
            if (IND)
                IND = IND(end);
                KSV(ii, 1:IND) = NaN;
                CHIABSV(ii, 1:IND) = NaN;
            end
        end
        FAMIN = -LAMV./(1-LAMV);

        figure('position', [0, 0, 600, 800]);
        subplot(2, 1, 1);hold
        surf(LAMV, FAV, CHIABSV, 'edgecolor', 'none', 'facecolor', 'interp')
        colormap copper
        view([0, 90]);grid off
        set(gca,'fontsize',18)

        xlabel('\lambda');ylabel('f_A')
        plot3(LAMV, 1-FAMIN, 1e5*ones(1, length(FAMIN)), 'k-', 'linewidth', 3)
        plot3(LAMV, FAMIN, 1e5*ones(1, length(FAMIN)), 'k-', 'linewidth', 3)
        ylim([0.1, 0.9]);box on
        colorbar
        caxis([0, 8])
        set(gca,'fontsize',18)

        subplot(2, 1, 2);hold
        surf(LAMV, FAV, KSV, 'edgecolor', 'none', 'facecolor', 'interp')
        colormap copper
        view([0, 90]);grid off
        set(gca,'fontsize',18)

        xlabel('\lambda');ylabel('f_A')
        plot3(LAMV, 1-FAMIN, 1e5*ones(1, length(FAMIN)), 'k-', 'linewidth', 3)
        plot3(LAMV, FAMIN, 1e5*ones(1, length(FAMIN)), 'k-', 'linewidth', 3)
        ylim([0.1, 0.9]);box on
        colorbar
        caxis([0, 4.5])
        set(gca,'fontsize',18)
    end

%     figure;
%     figure;hold
%     for ii = 1:51
%         COL = (ii-1)/(length(FAV)-1);
%         
%         FA = FAV(ii)
%         plot(LAMV, KSV(ii, :), ...
%             'o-', 'linewidth', 1.5, 'color', [COL 0 1-COL])
%     end
%     
%     figure;hold
%     for ii = 1:length(FAV)
%         COL = (ii-1)/(length(FAV)-1);
%         plot(LAMV, CHIABSV(ii,:)*FAV(ii)*(1-FAV(ii)), 'color', [COL 0 1-COL])
%     end
end