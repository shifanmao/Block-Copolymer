clear;

EPS = 0.01;
LAM = -0.75;
CHI = 0;

% simulation constants
FA=0.5;  % fraction of A blocks
M=8;        % number of blocks
G=5;
NM=G*EPS;  % number of Kuhn steps per monomer

R2=-0.5+0.5*exp(-2*NM)+NM;
k=logspace(-1,2,100)./sqrt(R2); % wavevectors

val = s2invwlc(M,NM,FA,LAM,k);
figure;
plot(k.*sqrt(R2),1./(-2*CHI+EPS*val),'k--','linewidth',3)