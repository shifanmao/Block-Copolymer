%this code tests the calculation of 3-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
function gam3 = gamma3(QM,ang,NM_WLC)

%Calculation parameters
ORDEig=10;
ORDL=10;
ResLayer=500;
ImagThreshold=1e-8;

%calculate gam3
gam3 = zeros(length(QM),length(ang));
for ii = 1:length(QM)
    for jj = 1:length(ang)
    fprintf('Calculating gam3 at QM = %.2f, ang = %.2f\n',QM(ii),ang(jj))

        Q1=QM(ii)*[1,0,0];
        Q2=transpose(rotz(ang(jj))*Q1');
        Q3=-Q1-Q2;

        gam3(ii,jj)=s3wlc(NM_WLC,Q1,Q2,Q3,ORDEig,ORDL,ResLayer);
        gam3(ii,jj)=gam3(ii,jj)/power(NM_WLC,3);

        if imag(gam3(ii,jj))<ImagThreshold
            gam3(ii,jj)=real(gam3(ii,jj));
        end
    end
end