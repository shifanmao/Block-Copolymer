%this code tests the calculation of 2-point correlation
clear;
addpath('../misc')

%Check ABs (if chkab==1, check combinations of SAB, SBA, etc.)
chkab=0;

%Chain structural information
NM=100;
N=100;

%Chain chemical information
FA=0.3;
LAM=-1;

% parameters for WLC calculations
ORDEig=10;
NumLayer=500;

%wavevector
k=logspace(-2,6,50)'/sqrt(r2(NM));

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% Gaussian Chain %%%%
%calculate s2
g2=zeros(length(k),2,2);
for ii=1:length(k)
    ii
%     g2(ii,:,:) = s2gc(N,NM_GS,LAM,FA,k(ii)/(NM_GS*N));
    g2(ii,:,:) = s2wlc(N,NM,LAM,FA,k(ii)/(NM*N));
%     g2(ii,:,:) = s2rigid(NM_GS,FA,k(ii)/NM_GS);
end
g2 = g2./power(NM*N,2);

leg=[];
%make plots
plot(k*sqrt(r2(NM)),g2(:,1,1),'b-',k*sqrt(r2(NM)),g2(:,2,2),'b--','linewidth',2);
leg=[leg {'Gaussian S_{AA}'}];
leg=[leg {'Gaussian S_{AB}'}];

legend(leg);
set(gca,'xscale','log');set(gca,'yscale','linear');ylim([0,1.2]);
xlabel('Normalized Wavevector kL');ylabel('Normalized Structure Factor S/L^2')