clear;

CHIV = 0.;
FA = 0.5;
Nvec = logspace(0,5,20);
Nbarvec = logspace(0,5,20);

chi1 = zeros(length(Nvec),length(Nbarvec));
for ii = 1:length(Nvec)
    for jj = 1:length(Nbarvec)
        [ii,jj]
        chi1(ii,jj)=spinodalRG(Nvec(ii),Nbarvec(jj),FA,CHIV);
    end
end

figure;hold
for ii = 1:length(Nvec)
    col = (ii-1)/(length(Nvec)-1);
    plot(Nbarvec,chi1(ii,:)*Nvec(ii),'color',[col 0 1-col]);
end
plot([1e0,1e5],[10.495,10.495],'k--')
%plot(Nbarvec,10.495+41.022*power(Nbarvec,-1/3),'k-')