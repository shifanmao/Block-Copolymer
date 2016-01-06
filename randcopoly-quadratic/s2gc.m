function [valaa,valab,valba,valbb]=s2gc(N,NM,FA,LAM,k)

% Calculate the Fourier transform of the Green function
% for the Gaussian chain chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Reset N to a row vector if entered as a column

if iscolumn(N)==1
    N=transpose(N);
end

valaa=zeros(length(k),length(N));
valab=zeros(length(k),length(N));
valba=zeros(length(k),length(N));
valbb=zeros(length(k),length(N));

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

for j=1:length(k)

    R=-k(j)^2/6;
    Z0=exp(R*NM);
    Z1=Z0*LAM;
    
    valeq=2*N*(Z0/R^2-1/R^2-NM/R);
    valne1=4/R^2*Z0*(Z0.^N-N*Z0+N-1)/(1-Z0)^2*(cosh(R*NM)-1);
    valne2=4/R^2*Z1*(Z1.^N-N*Z1+N-1)/(1-Z1)^2*(cosh(R*NM)-1);
    
    valeq(isnan(valeq))=0;
    valne1(isnan(valne1))=0;
    valne2(isnan(valne2))=0;
    
    valeq(isinf(valeq))=0;
    valne1(isinf(valne1))=0;
    valne2(isinf(valne2))=0;
    
    valaa(j,:)=valaa(j,:)+(FA*valeq+FA^2*valne1+FA*(1-FA)*valne2);
    valab(j,:)=valab(j,:)+(FA*(1-FA)*valne1-FA*(1-FA)*valne2);
    valbb(j,:)=valbb(j,:)+((1-FA)*valeq+(1-FA)^2*valne1+FA*(1-FA)*valne2);
    
    valba(j,:)=valab(j,:);
    
end