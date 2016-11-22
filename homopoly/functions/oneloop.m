%  calculate one-loop density-density correlation
clear;
%close all

u0 = 1;  % interaction parameter
c = 1;  % chain concentration

N=1e4;
k=logspace(0,5,100)';

gam2 = gamma2(k,N);
loop0 = power(1/u0+c*N^2*gam2,-1);

% perform integral
loop1Int = k.^2*power(gam2./(u0.^-1+c*gam2),2);


% % test zero wavevector limit of s3, i.e. s3(0,q) -> s2(q)
% % calculation parameters
% ORDEig=10;
% ORDL=10;
% ResLayer=500;
% ImagThreshold=1e-8;

% QM = logspace(-2,2,50)';
% s2 = zeros(length(QM),1);
% s3 = zeros(length(QM),1);
% ang = pi;
% for ii = 1:length(QM)
%     ii
%     Q1=QM(ii)*[1,0,0];
%     Q2=transpose(rotz(ang)*Q1');
%     Q3=[0,0,0];
%     s2(ii) = s2wlc(N,QM(ii),ORDEig,ResLayer)/power(N,2);
%     s3(ii) = s3wlc(N,Q1,Q2,Q3,ORDEig,ORDL,ResLayer)/power(N,3);
% end

% s = power(1./gam2+c*u*N^2,-1);
% 
% 
% figure;
% loglog(k,s)
% xlabel('k2l_p');ylabel('S(k)')