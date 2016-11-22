clear;
%close all

PLOTON = 1;
SAVEON = 1;

addpath('../functions/')
addpath('../chainstats/')
addpath(genpath('../chainstats/'))
addpath('../misc/')

% simulation constants
CHI=0;
FAV=0.5;  % fraction of A blocks
NV = 32*[5,2,1,.5,.2,.1,.05];

for N = NV;
    for FA = FAV
        fprintf('Calculating at N = %.2f, FA = %.2f\n', N, FA)

        % save destination
        filename = sprintf('sdata/SDIB_THRY_N%.2fFA%.2f',N,FA);

        % wavevectors
        %R2=-0.5+0.5*exp(-2*N)+N;
        k=logspace(-1,2,100); % wavevectors
        
        % calculate spinodal
        [chis,ks,~]=spinodal(N,FA);
        CHIS = chis*N;

        % calculate structure factor
        val = gamma2(N,FA,k,CHI);

        % save to result
        spin = [CHIS,ks];
        result = [k',val];  

        if (PLOTON)
            figure;
            plot(result(:,1),1./result(:,2),'--','linewidth',3)
            set(gca,'xscale','log');set(gca,'yscale','log');
        end

        if (SAVEON)
            dlmwrite(filename,[spin;result],'precision','%.6f');
        end
    end
end