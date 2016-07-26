function s2=s2wlc(NM,k,ORDEig,ResLayer)
% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimension
%
% Andrew Spakowitz (4/14/15)
s2=0;

% calculate the eigenvalues
R=Eigenvalues(k,ORDEig,1);
NR=ORDEig;

% get the residues for all roots
Residue=Residues(k,R(1:NR),NR,1,ResLayer,1);

for I=1:NR
    R(I)=R(I)*NM;
    valeq = R(I).^(-2).*((-1)+exp(1).^(R(I))+(-1).*R(I));
    s2=s2+2*Residue(I)*valeq*(NM^2);    
end

end