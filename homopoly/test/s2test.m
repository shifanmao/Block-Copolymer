%this code tests the calculation of 2-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
clear;
addpath(genpath('../chainstats/'))

%Calculation parameters
ORDEig=20;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
NM=1e2;

%wavevector
k=logspace(-2,2,100)';

%calculate s2
s2=zeros(length(k),1);
for ii=1:length(k)
    s2(ii) = s2wlc(NM,k(ii),ORDEig,ResLayer);            
    s2(ii) = s2(ii)/power(NM,2);
    
    if imag(s2(ii))<ImagThreshold
        s2(ii)=real(s2(ii));
    end
end

figure;hold;
plot(k,s2,'k-','linewidth',2);
xlim([1e0,1e3]);ylim([0,1]);box on

set(gca,'xscale','log');set(gca,'yscale','linear');
xlabel('Wavevector 2l_Pq');
ylabel('Pair Correlation $\left< \rho(\vec{q})\rho(-\vec{q}) \right>$','interpreter','latex')
xlim([1e-2,1e2])
set(gca,'fontsize',20)