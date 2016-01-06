%% Makes the plots of 'Impact of' paper
clear;close all;

PLOTON = 1;
EPSV = [0.01,0.10,1.00];
LAMV = [-0.75,-0.50,-0.25,0.00];

for EPS = EPSV
    for LAM = LAMV
        plotmcsim(EPS,LAM,PLOTON);
    end
end