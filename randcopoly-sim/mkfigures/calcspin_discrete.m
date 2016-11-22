% Some simple examples using functions to evaluate
% phase behavior of random copolymers of wormlike chains
% kmaxwlc.m and s2invwlc
%   -> use kmaxgc.m and s2invgc.m for Gaussian chains
%   -> use kmaxrr.m and s2invrr.m for perfectly rigid rods
%      (expect rigid rod functions to be relatively slow with
%       large number of monomers)
addpath('../functions')
close all
clear

% Simple Example 1: find spinodal, critical wavelength, and peak sharpness
NV=logspace(0,3,20);  % total of 100 monomers
NMV=logspace(-2,2,5); % each monomer has 0.1 Kuhn steps
FA=0.5; % equal chemical composition
LAM=-0.25;  % statistically anti-correlated random copolymer

figure(1);hold;set(gca,'fontsize',20)
figure(2);hold;set(gca,'fontsize',20)

cnt = 1;
for NM = NMV
    KS = [];
    CHIS = [];
    for N = NV
        fprintf('Calculating NM = 1e^%.2f, N = 1e%.2f \n',log10(NM),log10(N))
        [ks,sval,D2GAM2]=kmaxwlc(N,NM,FA,LAM);
        KS = [KS,ks];
        CHIS=[CHIS,0.5*sval];  % spinodal
    end

    fprintf('Results :\nSpinodal chivN_M= %.2f\n',CHIS*NM)
    fprintf('Critical wavemode q* = %.2f\n',KS)
    fprintf('Second der. of structure factor at q* = %.2f\n',D2GAM2)


    figure(1);plot(NV,CHIS*NM)
    set(gca,'xscale','log')

    figure(2);plot(NV,KS*sqrt(r2wlc(NM)))
    set(gca,'xscale','log')
end

figure(1);xlabel('N');ylabel('\chiN_M');box on
legend(num2str(NMV'))
filename = sprintf('chis_LAM%.2f',LAM);
saveas(gcf,strcat(filename,'.png'),'png')

figure(2);xlabel('N');ylabel('q^*R_M');box on
legend(num2str(NMV'))
filename = sprintf('ks_LAM%.2f',LAM);
saveas(gcf,strcat(filename,'.png'),'png')