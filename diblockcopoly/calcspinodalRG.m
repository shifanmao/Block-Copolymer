clear;

% parameters
FA = 0.5;
N = 1e5;
NMV = 1e3;
NbarV = logspace(0,5,10);

% results to return
chi1 = zeros(length(NMV),length(NbarV));

iNM=1;iNbar=1;
for NM = NMV
    for Nbar = NbarV
        Nbar
        chi1(iNM,iNbar)=spinodalRG(N,NM,Nbar,FA);
        iNbar = iNbar+1;
    end
    iNM = iNM+1;
end

figure;hold
plot([1e0,1e5],[10.495,10.495],'k--')
plot(NbarV,10.495+41.022*power(NbarV,-1/3),'k-')
plot(NbarV,chi1*N,'k+');
set(gca,'xscale','log')