function val=gamq2gc(M,NM,LAM,k)
%% Calculate the Fourier transform of the Green function
% for the Gaussian chain chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Inputs:
% M :: number of monomers
% NM :: number of Kuhn steps per monomer
% FA :: fraction of A
% LAM :: chemical correlation
% k :: wavevector normalized by 2l_p

%% Reset M to a row vector if entered as a column
%if iscolumn(M)==1
%    M=transpose(M);
%end

s=zeros(length(k),length(M));
for j=1:length(k)
    if k(j)*sqrt(r2(NM))<1e-2
        % zero wavemode limit
        sum = LAM*(-1-M*(LAM-1)+LAM^M)/(LAM-1)^2;
        s(j,:) = (M+2*sum)*NM^2;
    else
        R=-k(j)^2/6;
        Z0=exp(R*NM);
        Z1=Z0*LAM;

        valeq=2*M*(Z0/R^2-1/R^2-NM/R);
        valne=4/R^2*Z1*(Z1.^M-M*Z1+M-1)/(1-Z1)^2*(cosh(R*NM)-1);

        valeq(isnan(valeq))=0;
        valne(isnan(valne))=0;

        valeq(isinf(valeq))=0;
        valne(isinf(valne))=0;

        s(j,:)=s(j,:)+(valeq+valne);
    end
end
val=M./real(s);