function [val]=s2invrr(N,NM,FA,LAM,k)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Calculate the s matrix

[SAA,SAB,SBA,SBB]=s2rr(N,NM,FA,LAM,k);

DET=SAA.*SBB-SAB.*SBA;

val=real(N*NM*(SAA+SBB+2*SAB)./DET);
