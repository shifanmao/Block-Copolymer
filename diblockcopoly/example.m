% make a mean-field phase diagram
N=1e6;
FAV=linspace(0.1,0.499,41);
plotphase(N,FAV)

% make a phase diagram with density fluctuations
N=1e9;
Nbar=1e6;
FAV=linspace(0.3,0.499,41);
plotphaseRG(N,Nbar,FAV);