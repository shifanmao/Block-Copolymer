% Figure 2: Triangular phases at 8different 
clear;

CALCON = 1;
filename1 = '../data/figureS4ATriMesh';
filename2 = '../data/figureS4BTriMesh';
filename3 = '../data/figureS4CTriMesh';
filename4 = '../data/figureS4DTriMesh';

filename5 = '../data/figureS4ETriMesh';
filename6 = '../data/figureS4FTriMesh';
filename7 = '../data/figureS4GTriMesh';
filename8 = '../data/figureS4HTriMesh';

filenames = {filename1, filename2, filename3, filename4,...
             filename5, filename6, filename7, filename8};

if CALCON
    N = 100;
    cnt = 1;
%     for LAM = [-0.75, 0]
    for LAM = -0.75
        for NM = [1e-2, 1e2]
            for CHI = [0, 0.8]
                filename = filenames{cnt};
                
                [FAVV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, NM, LAM, CHI, 1, 0);
                save(filename, 'FAVV', 'PHIPV', 'EIGV', 'EIG', 'KSV', 'NM');
                cnt = cnt+1;
            end
        end
    end
else
    names = {};
    cnt = 1;
    for LAM = [-0.75, 0]
        for NM = [1e-2, 1e2]
            for CHI = [0, 0.8]
                name = sprintf('NM=%d, CHI=%.2f, LAM=%.2f', NM, CHI, LAM);
                names{cnt} = name;
                cnt = cnt+1;
            end
        end
    end
    
    for ii = 1:8
        filename = filenames{ii};
        load(filename);

        if ii <=4
            LAM = -0.75;
        else
            LAM = 0;
        end
        % put limits to FA
        FAMIN = -LAM/(1-LAM);
        FAMAX = 1/(1-LAM);
        for jj = 1:length(PHIPV)
            INDMIN = find(FAVV(:, jj) < FAMIN);
            INDMAX = find(FAVV(:, jj) > FAMAX);

            if INDMIN
                INDMIN = INDMIN(end);
                EIG(1:INDMIN, jj) = NaN;
            end
            if INDMAX
                INDMAX = INDMAX(1);
                EIG(INDMAX:end, jj) = NaN;
            end
        end

        figure;plotphase_wsolvent(FAVV, PHIPV, EIGV, EIG, KSV, NM)
        title(names{ii})
    end
end

% N = 1e3;
% NM = 1e2;
% FA = 0.5;
% PHIP = 1;
% LAM = 0;
% 
% % calculate spinodal CHIABS
% [CHIABS, KS, EIGS, EIGVS] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
% 
% RM = sqrt(r2(NM));
% KV = logspace(-1, 2, 100)/RM;
% 
% CHIABV = [0,0.2,0.4,0.6,0.8];
% figure;hold
% for ii = 1:length(CHIABV)
%     COL = (ii-1)/(length(CHIABV)-1)
% 
%     CHIAB = CHIABV(ii)*CHIABS;
%     CHIBA = CHIAB;
%     CHI = [0, CHIAB; CHIBA, 0];
% 
%     [EIG,EIGV]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,PHIP);
% 
%     plot(KV*RM, 1./EIG(:,1)/NM, '-','color', [COL 0 1-COL], 'linewidth', 1.5)
%     plot(KV*RM, 1./EIG(:,2)/NM, '--','color', [COL 0 1-COL], 'linewidth', 1.5)
% end
% 
% box on
% xlabel('Rq');
% %         ylabel('<\psi(-q)\psi(q)>')
% ylabel('$\langle \bf{\psi}^{\bf{\top}}(\vec{q}) \bf{\psi}(-\vec{q}) \rangle/Nv$',...
%     'Interpreter', 'Latex')
% xlim([min(KV), max(KV)]*RM)
% set(gca,'fontsize',18)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% %         title(strcat('N=', sprintf('%d', N), ', f_A=', sprintf('%.2f', FA)))
% 
% if N==1e3
%     legend('\chi_{AB}=0', '\chi_{AB}=0.2\chi_{AB}^{*}', '\chi_{AB}=0.4\chi_{AB}^{*}',...
%    '\chi_{AB}=0.6\chi_{AB}^{*}', '\chi_{AB}=0.8\chi_{AB}^{*}', 'location', 'northwest')
% else
%     legend('\chi_{AB}=0', '\chi_{AB}=0.2\chi_{AB}^{*}', '\chi_{AB}=0.4\chi_{AB}^{*}',...
%    '\chi_{AB}=0.6\chi_{AB}^{*}', '\chi_{AB}=0.8\chi_{AB}^{*}', 'location', 'southeast')
% end
% %         end
% %     end
% % end