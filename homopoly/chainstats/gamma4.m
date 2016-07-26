%this code tests the calculation of 4-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
function gam4 = gamma4(QM,ang,NM_WLC)
addpath('../chainstats')
addpath('../misc')
addpath('../eigcalc')

%Calculation parameters
ORDEig=8;
ORDL=4;
ResLayer=500;
ImagThreshold=1e-8;

%wavevector and structure factor
gam4 = zeros(length(QM),length(ang));
for ii = 1:length(QM)
    for jj = 1:length(ang)
        fprintf('Calculating Vd at QM = %.2f, ang = %.2f\n',QM(ii),ang(jj))
        Q1=QM(ii)*[1,0,0];
        Q2=transpose(rotz(ang(jj))*Q1');
        Q3=-Q1;
        Q4=-Q2;
        
        s4 = s4wlc(NM_WLC,Q1,Q2,Q3,Q4,ORDEig,ORDL,ResLayer);
        s4 = s4/power(NM_WLC,4);
        s2 = s2wlc(NM_WLC,QM(ii),ORDEig,ResLayer);            
        s2 = s2/power(NM_WLC,2);
    
        gam4(ii,jj) = s4-s2.^2;
        
        if imag(gam4(ii,jj))<ImagThreshold
            gam4(ii,jj)=real(gam4(ii,jj));
        end
    end
end