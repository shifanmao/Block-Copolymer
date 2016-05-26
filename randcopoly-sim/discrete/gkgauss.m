function [val]=gkgauss(N,G,LAM,FA,A1,A2,K)

% calculate the roots or eigenvalues of the Schrodinger equation

% k is a vector of all frequencies, for each k, get the roots

val=zeros(length(K),1);

NR=1;
Residue=1;

for j=1:length(K)
    
    val(j)=0;
    for I=1:NR
        LAMI=exp(-K(j)^2*G/6);
        PA1A2=pa1a2(A1,A2,FA,LAM,1,1);
        COM=PA1A2*N*Residue(I)*(G+2*LAMI*(LAMI^G-G*LAMI+G-1)/(1-LAMI)^2);
        COM(isnan(COM))=0;
        val(j)=val(j)+COM;
        
        for I2=2:N
            for I1=1:(I2-1)
                PA1A2=pa1a2(A1,A2,FA,LAM,I1,I2);
                COM=2*PA1A2*Residue(I)*(LAMI^(1-G)*(1-LAMI^G)^2/(1-LAMI)^2)*LAMI^(abs(I2-I1)*G);
                COM(isnan(COM))=0;
                val(j)=val(j)+COM;
            end
        end
    end
    
    if abs(imag(val(j)))<1e-14
        val(j)=real(val(j));
    end
end
