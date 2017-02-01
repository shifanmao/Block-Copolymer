function s2=s2wlc(N,NM,LAM,FA,k)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

FB=1-FA;
F=[FA,FB];
s2=zeros(2,2);

% if sum(power(k,2)) < 1e-8
%     s2(1,1)=
% end

% calculate the eigenvalues
ORDEig=20;  % maximum number of eigenvalues
ResLayer=500;  % number of residual layers
R=Eigenvalues(k,ORDEig,1);
NR=ORDEig;

% get the residues for all roots of each k(j)
Residue=Residues(k,R(1:NR),NR,1,ResLayer,1);
[j0,dj0]=CalRes0(k,ResLayer);
valeq=(NM/j0-dj0/j0^2);
s2(1,1)=s2(1,1)+2*N*F(1)*valeq;
s2(2,2)=s2(2,2)+2*N*F(2)*valeq;

for I=1:NR
    Z0=exp(R(I)*NM);
    Z1=Z0;
    ZL=Z0*LAM;
    valeq=Z0/R(I)^2;
    valne1=(2*Z1./N).*(Z1.^N-N*Z1+N-1)/(1-Z1)^2*(cosh(R(I)*NM)-1)/R(I)^2;
    valneL=(2*ZL./N).*(ZL.^N-N*ZL+N-1)/(1-ZL)^2*(cosh(R(I)*NM)-1)/R(I)^2;

    valeq(isnan(valeq))=0;
    valne1(isnan(valne1))=0;
    valneL(isnan(valneL))=0;

    valeq(isinf(valeq))=0;
    valne1(isinf(valne1))=0;
    valneL(isinf(valneL))=0;

    s2(1,1)=s2(1,1)+2*N*F(1)*Residue(I)*valeq;
    s2(2,2)=s2(2,2)+2*N*F(2)*Residue(I)*valeq;
    
    s2(1,1)=s2(1,1)+2*N*Residue(I)*(F(1)*F(1)*valne1 + F(1)*F(2)*valneL);
    s2(2,2)=s2(2,2)+2*N*Residue(I)*(F(2)*F(2)*valne1 + F(1)*F(2)*valneL);
    s2(1,2)=s2(1,2)+2*N*Residue(I)*(F(1)*F(2)*valne1 - F(1)*F(2)*valneL);
end
s2(2,1)=s2(1,2);

s2(imag(s2)<1e-2) = real(s2);
    
end

function [j0,dj0]=CalRes0(k,ResLayer)
d = 3;

if k<=1
    
    % residual using continued fraction for small k
    p=0;
    W=p+(ResLayer+d-2)*ResLayer;
    Wprime=1;
    for L=ResLayer:-1:1
        AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
        Wprime=1-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm+AL^2/W;
    end
    dj0=Wprime;
    j0=W;
else
    
    % residual using continued fraction for large k
    p=0;
    W=(p+(ResLayer+d-2)*ResLayer)/k;
    Wprime=1/k;
    for L=ResLayer:-1:1
        AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
        Wprime=1/k-AL^2*Wprime/W^2;
        PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
        W=PLm/k+AL^2/W;
    end
    dj0=k*Wprime;
    j0=k*W;
end

end