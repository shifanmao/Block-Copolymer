function randphase_func(EPS)

N=20;
G=5;
FA=0.5;

LAM0=-1;
LAMF=1;
NLAM=2001;
LAMV=transpose(linspace(LAM0,LAMF,NLAM));
%LAMV=-0.75:0.25:0.75;NLAM=length(LAMV);

EPS0=1e-2;
EPSF=1e3;
NEPS=51;
EPSV=transpose(logspace(log10(EPS0),log10(EPSF),NEPS))/G;
%EPSV=transpose(linspace(EPS0,EPSF,NEPS));

ORD=30;         % (Number of poles-1) in Laplace inversion
ResLayer=100;   % Number of layers used in the continued fraction form
ORDmax=30;

R2=zeros(NLAM,NEPS);
CHI=zeros(NLAM,NEPS);
KS=zeros(NLAM,NEPS);

test=[];
%for IE=1:NEPS
%    EPS=EPSV(IE);
    NG=EPS*G;
    for IL=1:NLAM
        R2(IL)=r2wlc(NG);
        
        LAM=LAMV(IL);
        
        [KMAX,SMAX]=kmax(N,G,LAM,FA,EPS,ORDmax,ORD,ResLayer);
        
        KS(IL)=KMAX;
        CHI(IL)=SMAX/2;
        
        %[IL,KMAX,SMAX/2]
        [NG,LAM,SMAX/2*G]
        test=[test;NG,LAM,KMAX,sqrt(r2wlc(NG)),SMAX/2*G];
        filename=sprintf('../Sdata/spindata');
        dlmwrite(filename,test,'delimiter','\t','precision',6)
    end
    %save('data/phase','R2','CHI','KS','LAMV','EPSV','FA','N','G')
%end
