%this code tests the calculation of 3-point correlation
clear;
addpath('misc')

%Chain structural information
NM_GS=1;
N=100;

%Chain chemical information
FA=0.1;LAM=-0.75;

% parameters for WLC calculations
ORDEig=5;
ORDL=4;
NumLayer=500;

%wavevector and structure factor
QM=logspace(0,3,20)';
% QM = 5*N*NM_GS;
Q1=zeros(length(QM),1);
Q2=zeros(length(QM),1);
Q3=zeros(length(QM),1);
ang=2*pi/3;
for ii=1:length(QM)
    Q1(ii,1:3)=QM(ii)*[1,0,0]';
    Q2(ii,1:3)=transpose(rotz(ang)*Q1(ii,1:3)');
    Q3(ii,1:3)=-Q1(ii,1:3)-Q2(ii,1:3);
end

%begin making plots
figure;hold;set(gca,'fontsize',15);leg=[];

%%%% Gaussian Chain %%%%
%calculate s3
g3=zeros(length(QM),2,2,2);
for ii=1:length(QM)
    ii
%     g3(ii,:,:,:)=s3gaussian(N,NM_GS,LAM,FA,...
%                                      Q1(ii,1:3)/(N*NM_GS),...
%                                      Q2(ii,1:3)/(N*NM_GS),...
%                                      Q3(ii,1:3)/(N*NM_GS));
    g3(ii,:,:,:)=s3wlc(N,NM_GS,LAM,FA,...
                                Q1(ii,1:3)/(N*NM_GS),...
                                Q2(ii,1:3)/(N*NM_GS),...
                                Q3(ii,1:3)/(N*NM_GS),ORDEig,ORDL,NumLayer);
%     g3(ii,:,:,:)=s3rigid(NM_GS,FA,Q1(ii,1:3)/NM_GS,...
%                                 Q2(ii,1:3)/NM_GS,...
%                                 Q3(ii,1:3)/NM_GS);
end
g3=g3/power(NM_GS*N,3);    

%make plots
plot(QM,g3(:,1,1,1),'k-',...
     QM,g3(:,1,2,1),'b-',...
     QM,g3(:,1,1,2),'r-','linewidth',2);
leg=[leg {['Gaussian S_{AAA}']}];
leg=[leg {['Gaussian S_{AAB}']}];

set(gca,'xscale','log');set(gca,'yscale','linear');