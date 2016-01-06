function [dum]=r2calcfun(testnum,chemcp)

  %%%%%%%%% PARAMS %%%%%%%%%%%%
  boxl=20;dk=0.2;Ree=2.0;
  FA=0.5;
  % EPS=0.01;LAM=-0.75;
  G=5;N=8;
  %%%%%%%%% END OF PARAMS %%%%%%

  %%%%%%%%% INPUT %%%%%%%%%%
  chilistname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/',testnum,chemcp);
  SPINODAL=load('/tower12/home/shifan/polymem-wlc-pt/scalcbatch/analytical/chivals');
  CHIV=load([chilistname,'chilist']);
  CHIV=CHIV(1:5:46);
  NCHI=length(CHIV);

  snap0=1;snapf=101;skip=1;
  snapv=snap0:skip:snapf;
  nsnap=length(snapv);
  %%%%%%%%% INPUT END %%%%%%%%%%

  %%%% start calculations %%%%
  rgvec = zeros(NCHI,1);  % radius of gyration
  r2vec = zeros(NCHI,1);  % end-to-end separation of a monomer
  for cnum=1:NCHI
    CHI=CHIV(cnum);
    col=(cnum-1)/(NCHI-1);

    test=[];
    rg = 0;
    r2 = 0;
    for kk=snapv
      disp(['chi=',num2str(CHI),' // snap#', num2str(kk)])
      cmat=rankpt(testnum,chemcp,kk);
      ind=find(abs(cmat(:,2)-CHI)<1e-2);
      index=cmat(ind,1);
      indexname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/rand-wlc-%d/',testnum,chemcp,index);

      r=load([indexname,'data/r',num2str(kk)]);
      [a,b] = rcalc(r,N,G);
      rg = rg + a;
      r2 = r2 + b;
    end
    rgvec(cnum) = rg./nsnap;
    r2vec(cnum) = r2./nsnap;
  end

  test=[CHIV,rgvec,r2vec];
  filename=sprintf('rdata/RMC_SIM%dCHEM%d',testnum,chemcp);
  dlmwrite(filename,test,'delimiter','\t','precision',3)
end

function [rg,r2] = rcalc(r,N,G)
  rg = 0;
  r2 = 0;
  cnt = 1;
  for ii = 1:N*G:length(r)
    rg = rg + rgcalc(r(ii:ii+N*G-1,1:3));
    r2 = r2 + r2calc(r(ii:ii+N*G-1,1:3),G);
    cnt = cnt+1;
  end
  rg = rg/cnt;
  r2 = r2/cnt;
end

function rg = rgcalc(r)
  % calculates radius of gyration
  rmean = mean(r,1);
  N = length(r);
  rg = 0;
  for ii = 1:N
    rg = rg+(1/N)*norm(r(ii,:)-rmean)^2;
  end
  rg = sqrt(rg);
end

function r2 = r2calc(r,G)
  % calculates end-to-end separation of monomer
  r2 = 0;
  cnt = 1;
  for ii = 1:G*2:length(r)
    % r2 = r2 + norm(r(ii+G-1)-r(ii));
    rc1 = mean(r(ii:ii+G-1,:));
    rc2 = mean(r(ii+G:ii+2*G-1,:));
    r2 = r2 + norm(rc1-rc2);
  end
  r2 = r2/cnt;
  r2 = sqrt(r2);
end
