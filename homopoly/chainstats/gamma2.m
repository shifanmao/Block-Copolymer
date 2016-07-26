%this code tests the calculation of 2-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
function gam2 = gamma2(k,NM_WLC)
addpath('../chainstats/')
addpath('../eigcalc/')

%Calculation parameters
ORDEig=20;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
% NM_WLC=1e4;
% Rg = sqrt(NM_WLC/6);

%wavevector
% k=logspace(0,5,100)'/Rg;

%calculate s2
gam2=zeros(length(k),1);
for ii=1:length(k)
    s2 = s2wlc(NM_WLC,k(ii),ORDEig,ResLayer)/power(NM_WLC,2);
    
    if imag(s2)<ImagThreshold
        gam2(ii)=real(s2);
    end
end
