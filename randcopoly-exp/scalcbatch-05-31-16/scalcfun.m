function [dum]=scalcfun(testnum,chemcp,lk,vfac)

simdir = 'polymem-exp-05-31-16/rand-pt-exp-05-31-16'
calcdir = 'polymem-exp-05-31-16/scalcbatch-05-31-16'

%%%%%%%%% PARAMS %%%%%%%%%%%%
vfac = 2.; %scaling parameter
lfac = vfac^(1/3);
boxl = 20/lfac;

%readjust box size according to no. of polymers
N=8;G=5*vfac;v = 0.1/vfac;
NP=round(boxl^3/(N*G*v));
boxl=(v*N*G*NP)^(1/3);

Ree=2.0;lksample=20;
%%%%%%%%% END OF PARAMS %%%%%%

%%%%%%%%% INPUT %%%%%%%%%%
chilistname=strcat('/tower12/home/shifan/',simdir,sprintf('-%d-%d/',testnum,chemcp));
SPINODAL=load('/tower12/home/shifan/polymem-wlc-pt/scalcbatch/analytical/chivals');

CHIV=load([chilistname,'chilist'])
%CHIV=CHIV([1:2:11,21,31,41]);

NCHI=length(CHIV);

snap0=1;snapf=101;skip=1;
snapv=snap0:skip:snapf;
nsnap=length(snapv);
%%%%%%%%% INPUT END %%%%%%%%%%

%%%% generate vectors on unit sphere %%%%
[kb,numk]=basisgen(lk,lksample,boxl);nk=length(kb);

%%%% start calculations %%%%
for cnum=1:NCHI
CHI=CHIV(cnum);
col=(cnum-1)/(NCHI-1);

test=[];knum=1;
k=[];S=[];
Stot=zeros(nk,1);

for kk=snapv
disp(['chi=',num2str(CHI),' // snap#', num2str(kk)])
cmat=rankpt(testnum,chemcp,kk);
ind=find(abs(cmat(:,2)-CHI)<1e-2);
index=cmat(ind,1);
%indexname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/rand-wlc-%d/',testnum,chemcp,index);
indexname=strcat('/tower12/home/shifan/',simdir,sprintf('-%d-%d/rand-wlc-%d/',testnum,chemcp,index));

    r=load([indexname,'data/r',num2str(kk)]);
    n=length(r);

    xr=r(:,1);yr=r(:,2);zr=r(:,3);
    xr=xr-boxl*floor(xr./boxl);
    yr=yr-boxl*floor(yr./boxl);
    zr=zr-boxl*floor(zr./boxl);

    id=r(:,4);
    f=sum(r(:,4)/n);

    % transform chemical identity to \sigma_A-f_A
    t=2*(id-0.5);
    t=0.5*(t+1-2*f);

    for ii=1:nk
      xcom=kb(ii,1)*xr;
      ycom=kb(ii,2)*yr;
      zcom=kb(ii,3)*zr;

      Scos=(cos(xcom+ycom+zcom).*t);
      Ssin=(sin(xcom+ycom+zcom).*t);
      Stot(ii)=sum(Scos).^2+sum(Ssin).^2;
    end

    [k,S(knum,:)]=scalcavg(kb,Stot);
    knum=knum+1;
end

S=S./n;  %normalization
Sa=mean(S,1);

test=[k*Ree,transpose(Sa)];
filename=sprintf('../Sdata/SMC_SIM%dCHEM%dCHI%.8f',testnum,chemcp,CHI);
dlmwrite(filename,test,'delimiter','\t','precision',10)

end
