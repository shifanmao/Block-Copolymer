%this code tests the calculation of 4-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
clear;
addpath('misc')

%Calculation parameters
ORDEig=2;
ORDL=1;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
NM=1000;

%wavevector and structure factor
QM=logspace(-1,2,20)'/NM;
Q1=zeros(length(QM),3);
Q2=zeros(length(QM),3);
Q3=zeros(length(QM),3);
Q4=zeros(length(QM),3);
ang=pi/2;
for ii=1:length(QM)
    Q1(ii,1:3)=QM(ii)*[1,0,0];
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q1(ii,1:3);
    Q4(ii,1:3)=-Q2(ii,1:3);
end

%calculate s4
s4=zeros(length(QM),1);
for ii=1:length(QM)
    ii
    s4(ii)=s4wlc(NM,Q1(ii,1:3),Q2(ii,1:3),Q3(ii,1:3),Q4(ii,1:3),...
                        ORDEig,ORDL,ResLayer);
    s4(ii)=s4(ii)/power(NM,4);
end

figure;plot(QM*NM,s4,'b-','linewidth',2);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^4')