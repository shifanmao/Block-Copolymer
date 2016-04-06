function out = powl(m,x,N)
fact=[1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800];

if m==0
    out=x^N;
    return
end
if abs(x-1) < 0.01
    prodic=ones(1,12);
    for j=2:12
        prodic(j)=prodic(j-1)*(N-j+2);
    end
    out=sum(((x-1).^[m:11].*prodic(m+1:12)./fact(m+1:12)));
else
    prodic=ones(1,m);
    for j=2:m
        prodic(j)=prodic(j-1)*(N-j+2);
    end    
    out=x^N-sum(((x-1).^[0:m-1]).*prodic./fact(1:m));
end
end