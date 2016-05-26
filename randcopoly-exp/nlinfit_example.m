clear;

% Generate X vectors for both data sets
x1 = 0:0.1:10;
x2 = 0:1:10;

% Generate Y data with some noise
y1 = cos(2*pi*0.5*x1).*exp(-x1/5) + 0.05*randn(size(x1));
y2 = 0.5 + 2*exp(-(x2/5)) + 0.05*randn(size(x2));

% Define fitting functions and parameters, with identical
% exponential decay for both data sets
mdl1 = @(beta,x) cos(2*pi*beta(1)*x).*exp(-x/beta(2));
mdl2 = @(beta,x) beta(4) + beta(3)*exp(-(x/beta(2)));

% Prepare input for NLINMULTIFIT and perform fitting
x_cell = {x1, x2};
y_cell = {y1, y2};
mdl_cell = {mdl1, mdl2};
beta0 = [1, 1, 1, 1];
[beta,r,J,Sigma,mse,errorparam,robustw] = ...
            nlinmultifit(x_cell, y_cell, mdl_cell, beta0);

% Calculate model predictions and confidence intervals
[ypred1,delta1] = nlpredci(mdl1,x1,beta,r,'covar',Sigma);
[ypred2,delta2] = nlpredci(mdl2,x2,beta,r,'covar',Sigma);

% Calculate parameter confidence intervals
ci = nlparci(beta,r,'Jacobian',J);

% Plot results
figure;
hold all;
box on;
scatter(x1,y1);
scatter(x2,y2);
plot(x1,ypred1,'Color','blue');
plot(x1,ypred1+delta1,'Color','blue','LineStyle',':');
plot(x1,ypred1-delta1,'Color','blue','LineStyle',':');
plot(x2,ypred2,'Color',[0 0.5 0]);
plot(x2,ypred2+delta2,'Color',[0 0.5 0],'LineStyle',':');
plot(x2,ypred2-delta2,'Color',[0 0.5 0],'LineStyle',':');