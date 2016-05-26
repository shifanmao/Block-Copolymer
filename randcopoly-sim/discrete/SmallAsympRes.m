function Res=SmallAsympRes(K,n)

% calculate the residual using small k asymptot
n=n-1;
Res=1;
Wn=-n*(n+1);

for j=1:n
    Wjm1=-(j-1)*j;
    aj=j/sqrt(4*j^2-1);
    Res=Res*aj/(Wn-Wjm1);
end

Res=Res^2*K^(2*n)*(-1)^n;