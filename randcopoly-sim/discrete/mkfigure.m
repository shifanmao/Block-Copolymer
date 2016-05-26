clear;
close all;

%%%%%%%%%%%%%%% PLOT 1: Spinodal %%%%%%%%%%%%%%%%%%
data=load('data/N8G5_2');
EPSV=logspace(-2,0,13);
%EPSV=EPSV(1:3);
EPSV=fliplr(EPSV);
NEPS=length(EPSV);

figure;hold
set(gca,'fontsize',18);
for ii=1:NEPS
    col=(ii-1)/(NEPS-1);
    
    G=5;EPS=EPSV(ii);
    NG=EPS*G;
    
    ind1=find(abs(data(:,1)-EPSV(ii))<1e-3);
    %plot(data(ind1,2),2*pi./(data(ind1,3)./sqrt(r2wlc(NG))),'color',[col 0 1-col],'linewidth',1.2);
    plot(data(ind1,2),data(ind1,3),'color',[col 0 1-col],'linewidth',1.2);
end
xlim([-1,1]);ylim([0,5])
xlabel('\lambda');ylabel('Critial wavelength k^* R_b')

%%%%%%%%%%%%%%% PLOT 2: k modes %%%%%%%%%%%%%%%%%%
data=load('data/N8G5');
EPSV=logspace(-2,0,13);
%EPSV=EPSV(1:3);
EPSV=fliplr(EPSV);
NEPS=length(EPSV);

figure;hold
set(gca,'fontsize',18);
for ii=1:NEPS
    col=(ii-1)/(NEPS-1);
    ind1=find(abs(data(:,1)-EPSV(ii))<1e-3);
    plot(data(ind1,2),data(ind1,4),'color',[col 0 1-col],'linewidth',1.2);
end
xlim([-1,1]);ylim([0,6])
xlabel('\lambda');ylabel('Flory-Huggins Parameter 4f_A (1-f_A) \chi N_b')

%%%%%%%%%%scaling of D%%%%%%%%%
data=load('data/N8G5_2');
LAM0=-1;
LAMF=0;
LAMV=transpose(linspace(LAM0,LAMF,81));
LAMV=LAMV(1:10:51);
NLAM=length(LAMV);

figure;hold
set(gca,'fontsize',18);
for ii=1:NLAM
    col=(ii-1)/(NLAM-1);
    
    ind1=find(abs(data(:,2)-LAMV(ii))<1e-3);
    
    %G=5;EPS=data(ind1,2);
    %NG=EPS*G;
    %R2=r2wlc(NG)
    
    %2*pi/k vs. EPS*G
    plot(data(ind1,1)*G,2*pi./(data(ind1,3)./sqrt(r2wlc(data(ind1,1)*G))),...
                                    'color',[col 0 1-col],'linewidth',1.2);
end
xlabel('N_b');ylabel('D/2l_P')
set(gca,'xscale','log');set(gca,'yscale','log');
title('blue = anticorrelated, red = more correlated')