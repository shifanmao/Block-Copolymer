
x=0.9;
N=5.5;
tic

for count=1:1000
    prodic=ones(1,12);
    for j=2:12
        prodic(j)=prodic(j-1)*(N-j+2);
    end
    
end  
toc

prodic

for count=1:1000
    prodic=
    
end