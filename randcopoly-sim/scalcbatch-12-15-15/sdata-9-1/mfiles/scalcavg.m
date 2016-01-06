function [k,S]=scalcavg(kb,Stot)

tol=1e-3;
k=round(kb/tol)*tol;
k=unique(k(:,4));

numk=length(k);
S=zeros(numk,1);

for ii=1:numk
  index=find(abs(kb(:,4)-k(ii))<=tol);
  S(ii)=mean(Stot(index));
end
