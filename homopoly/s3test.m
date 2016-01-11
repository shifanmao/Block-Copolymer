%this code tests the calculation of 3-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
clear;

%Calculation parameters
ORDEig=10;
ORDL=10;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
NM_WLC=10;

%wavevector and structure factor
QM=logspace(-1,3,50)';
Q1=zeros(length(QM),1);
Q2=zeros(length(QM),1);
Q3=zeros(length(QM),1);
ang=pi;
for ii=1:length(QM)
    Q1(ii,1:3)=QM(ii)*[1,0,0];
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q1(ii,1:3)-Q2(ii,1:3);
end

%calculate s3
s3=zeros(length(QM),1);
for ii=1:length(QM)
    s3(ii)=s3wlc(NM_WLC,Q1(ii,1:3)/NM_WLC,...
                        Q2(ii,1:3)/NM_WLC,...
                        Q3(ii,1:3)/NM_WLC,ORDEig,ORDL,ResLayer);
    s3(ii)=s3(ii)/power(NM_WLC,3);
end

figure;plot(QM,s3,'b-','linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^3')