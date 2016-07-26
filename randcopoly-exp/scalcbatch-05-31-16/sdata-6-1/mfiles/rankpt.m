function cmat=rankpt(testnum,chemcp,runnum)

simdir = 'polymem-exp-05-31-16/rand-pt-exp-05-31-16'
%replica span
snum1=1;
snum2=22;
snum=snum2-snum1+1;
cmat=zeros(snum,2);

%snapshot span
dnum1=1;
dnum2=runnum;

dnum=dnum2-dnum1+1;
xsteps=dnum1:dnum2;
d=zeros(snum,dnum);G=5;

for ii=snum1:snum2
    indexname=strcat('/tower12/home/shifan/',simdir,sprintf('-%d-%d/rand-wlc-%d/',testnum,chemcp,ii-1));
    test=load([indexname,'data/out2']);
    cmat(ii,1)=ii-1;
    cmat(ii,2)=test(dnum2,2);
end
