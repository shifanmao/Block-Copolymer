%this code tests the calculation of 4-point correlation
%functions of rigid-rod, wormlike chain, and Gaussian chains
clear;
addpath(genpath('../chainstats'))

%Calculation parameters
ORDEig=10;
ORDL=10;
ResLayer=500;
ImagThreshold=1e-8;

%Number of Kuhn steps
NM=1;

%wavevector and structure factor
k=logspace(-2,3,21)';
Q1=zeros(length(k),3);
Q2=zeros(length(k),3);
Q3=zeros(length(k),3);
Q4=zeros(length(k),3);
ang=pi;
for ii=1:length(k)
    Q1(ii,1:3)=k(ii)*[1,0,0];
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q1(ii,1:3);
    Q4(ii,1:3)=-Q2(ii,1:3);
end

%calculate s4
s4=zeros(length(k),1);
for ii=1:length(k)
    ii
    s4(ii)=s4wlc(NM,Q1(ii,1:3),Q2(ii,1:3),Q3(ii,1:3),Q4(ii,1:3),...
                        ORDEig,ORDL,ResLayer);
    s4(ii)=s4(ii)/power(NM,4);
end

%calculate s2
s2=zeros(length(k),1);
for ii=1:length(k)
    s2(ii) = s2wlc(NM,k(ii),ORDEig,ResLayer);            
    s2(ii) = s2(ii)/power(NM,2);
    
    if imag(s2(ii))<ImagThreshold
        s2(ii)=real(s2(ii));
    end
end

% figure;plot(k,s4,'b-','linewidth',2);
figure; plot(k,s4,'kx',k,s2.^2,'k-')
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Wavevector 2l_Pq');
ylabel('Normalized Structure Factor S/L^4')
legend('s4','s2')