function s2=s2wlc(N,NM,FA,k)
%% function s2wlc :: Calculate the Fourier transform of the Green function
% for the worm-like chain in d-dimension
% Usage: s2=s2wlc(NM,FA,k)
% Inputs ::
%   N, Number of monomers
%   NM, Number of Kuhn steps per monomer
%   FA, Fraction of A monomers
%   k, wavevector
% Andrew Spakowitz (4/14/15)

addpath('eigcalc/')
% parameters for worm-like chain calculations
% d=3;  % number of dimensions
ORDmax=20;  % maximum number of eigenvalues
ResLayer=500;  % number of residual layers
FB=1-FA;

% calculate the eigenvalues
R=Eigenvalues(k,ORDmax,1);
NR=ORDmax;

% get the residues for all roots of each k(j)
Residue=Residues(k,R(1:NR),NR,1,ResLayer,1);

s2=zeros(2,2);
for I=1:NR
% Case 1 :: A1==A2
    R(I)=R(I)*NM;
    
    % on same monomer (S integrals)
    valeqA = R(I).^(-2).*(-1+exp(FA.*R(I))-FA.*R(I));
    valeqB = R(I).^(-2).*(-1+exp(FB.*R(I))-FB.*R(I));

    s2(1,1)=s2(1,1)+2*Residue(I)*valeqA*(N^2);
    s2(2,2)=s2(2,2)+2*Residue(I)*valeqB*(N^2);

% Case 2 :: A1~=A2
    % on same monomer (S integrals)
    valeqAB = (1+exp(R(I))-exp(FA.*R(I))-exp(R(I)-FA.*R(I))).*R(I).^(-2);
    % J1<J2
    s2(1,2)=s2(1,2)+Residue(I)*valeqAB*(N^2);
    
end
s2(2,1)=s2(1,2);

end