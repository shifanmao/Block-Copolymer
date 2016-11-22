clear;
%close all

PLOTON = 1;
SAVEON = 0;

addpath('../functions/')
addpath('../misc/')

% simulation constants
CHI=0;
M=100;     % number of blocks
G=5;  % number of discrete monomers

FAV=0.5;  % fraction of A blocks
EPSV = [0.01,0.10,1.00];
LAMV = [-1.00:0.25:0.75];

EPSV = 0.01;LAMV = -1;

% EPSV = [0.01,0.10,1.00,0.04];
% LAMV = 0;
% FAV = [0.16,0.23,0.50];

for EPS = EPSV
    for LAM = LAMV;
        for FA = FAV
            fprintf('Calculating at EPS = %.2f, LAM = %.2f, FA = %.2f\n', EPS, LAM, FA)
            NM=G*EPS;  % number of Kuhn steps per monomer

            % save destination
            filename = sprintf('sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);

            % wavevectors
            R2=-0.5+0.5*exp(-2*NM)+NM;
            k=logspace(-1,2,100)./sqrt(R2); % wavevectors

            % calculate spinodal
            [kval,sval]=kmaxwlc(M,NM,FA,LAM);
            CHIS=0.5*sval;

            % calculate structure factor
            val = s2invwlc(M,NM,FA,LAM,k);

            % save to result
            spinodal = [CHIS*NM,kval*sqrt(R2)];
            result = [k'.*sqrt(R2),val];

            if (PLOTON)
                figure;
                plot(result(:,1),1/EPS./result(:,2),'--','linewidth',3)
                set(gca,'xscale','log');set(gca,'yscale','log');
            end

            if (SAVEON)
                dlmwrite(filename,[spinodal;result],'precision','%.6f');
            end
        end
    end
end