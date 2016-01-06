function cmat=rankpt(testnum,chemcp,runnum)

%testnum=3;
%replica span
snum1=1;
snum2=46;
snum=snum2-snum1+1;
cmat=zeros(snum,2);

%snapshot span
dnum1=1;
dnum2=runnum;

dnum=dnum2-dnum1+1;
xsteps=dnum1:dnum2;
d=zeros(snum,dnum);G=5;

for ii=snum1:snum2
test=load(['/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-',num2str(testnum),'-',num2str(chemcp),'/rand-wlc-',num2str(ii-1),'/data/out2']);
    cmat(ii,1)=ii-1;
    cmat(ii,2)=test(dnum2,2);
end
