function s2inv=s2inverse(N,NM,LAM,FA,k)
%% calculates the inverse of pair correlation function in Fourier space
% usage s2inv=s2inverse(N,FA,k)
% Return:
%    val, inverse s2
% Parameters:
%    N, number of Kuhn steps per monomer
%    FA, fraction of A type monomers
%    k, magnitude of wavevector in unit of 1/contour length
%       k input can be a vector
s2inv=zeros(2,2);
MIN=1e-10;

% Calculate the s matrix
if NM>=1e4  % Gaussian chain limit
    s2=s2gc(N,NM,LAM,FA,k);
elseif NM<=1e-4  % Rigid rod limit
    s2=s2rr(N,NM,LAM,FA,k);
else
    s2=s2wlc(N,NM,LAM,FA,k);
end
DET=s2(1,1,:)*s2(2,2,:)-s2(1,2,:)*s2(2,1,:);

% if abs(k)<=MIN
%     s2inv=ones(2,2)./power(N*NM,2);
% else
    s2inv(1,1) = s2(2,2)./DET;
    s2inv(1,2) = -s2(1,2)./DET;
    s2inv(2,1) = -s2(2,1)./DET;
    s2inv(2,2) = s2(1,1)./DET;
end
