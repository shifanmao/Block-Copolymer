function f = saxsfitfunc(x)
global sf qf N rm FA NFIT

NQ = length(qf);
f = zeros(NQ*NFIT,1);

SCALE = x(1);
LAM = x(2);
NM = x(3);
R2 = r2wlc(NM);

for IT = 1:NFIT
    CHI = x(3+IT);
    
    fitfunc = SCALE./(-2*CHI+s2invwlc(N,NM,FA,LAM,qf*rm/sqrt(R2)));
    fitfunc = fitfunc./NM;
    %f(1+(IT-1)*NQ : IT*NQ) = fitfunc-sf(:,IT);
    f(1+(IT-1)*NQ : IT*NQ) = (fitfunc-sf(:,IT))./sf(:,IT);
end