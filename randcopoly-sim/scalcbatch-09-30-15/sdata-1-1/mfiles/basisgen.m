function [kbnew,nknew]=basisgen(N,Nsample,L)
%function [kb,nk]=basisgen(N,Nsample,L)
%Nsample=20;  %total number of bins
%L=10;  %size of simulation

kb=[];

for n1=0:Nsample
  for n2=0:Nsample
    for n3=0:Nsample
      if (n1+n2+n3)>0
        k1=(2*pi/L)*n1;
        k2=(2*pi/L)*n2;
        k3=(2*pi/L)*n3;
        km=sqrt(k1*k1+k2*k2+k3*k3);
        kb=[kb;k1,k2,k3,km];
      end
    end
  end
end

nk=length(unique(kb(:,4)));

%only select some k values at high k
kmag=unique(kb(:,4));
ind =unique(round(logspace(0,log10(nk),N)));
whos ind

ind2=[];kbnew=[];
for ii=1:length(ind)
  ind2=find(kb(:,4)==kmag(ind(ii)));
  kbnew=[kbnew;kb(ind2,:)];
end
nknew=length(unique(kbnew(:,4)));
