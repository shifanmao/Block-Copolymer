%this code tests the calculation of 2-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
clear;

%Calculation parameters
ORDEig=20;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
NM_WLC=100;

%wavevector
k=logspace(-1,3,100)';

%calculate s2
s2=zeros(length(k),1);
for ii=1:length(k)
    s2(ii) = s2wlc(NM_WLC,k(ii)/NM_WLC,ORDEig,ResLayer);            
    s2(ii) = s2(ii)/power(NM_WLC,2);
    
    if imag(s2(ii))<ImagThreshold
        s2(ii)=real(s2(ii));
    end
end

figure;plot(k,s2,'linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^2')