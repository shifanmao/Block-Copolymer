function [val]=gk(N,G,LAM,FA,EPS,A1,A2,K,ORDmax,ORD,ResLayer)

% calculate the roots or eigenvalues of the Schrodinger equation

% k is a vector of all frequencies, for each k, get the roots

val=zeros(length(K),1);

for j=1:length(K)
    
    % calculate the eigenvalues
    KJ=abs(K(j));
    
    if KJ ==0
        val(j)=0;
        PA1A2=pa1a2(A1,A2,FA,LAM,1,1);
        COM=PA1A2*N*G^2;
        COM(isnan(COM))=0;
        val(j)=val(j)+COM;
        
        for I1=1:(N-1)
            NUM=N-I1;
            PA1A2=pa1a2(A1,A2,FA,LAM,I1,0);
            COM=NUM*2*PA1A2*G^2;
            COM(isnan(COM))=0;
            val(j)=val(j)+COM;
        end
        
    else
        
        R=MatRoots(KJ,ORD);
        NR=ORDmax;
        
        % get the residues for all roots of each k(j)
        Residue=NewCalRes(ResLayer,R(1:NR),KJ);
        
        val(j)=0;
        for I=1:NR
            LAMI=exp(R(I)*EPS);
            PA1A2=pa1a2(A1,A2,FA,LAM,1,1);
            COM=PA1A2*N*Residue(I)*(G+2*LAMI*(LAMI^G-G*LAMI+G-1)/(1-LAMI)^2);
            COM(isnan(COM))=0;
            val(j)=val(j)+COM;
            
            for I1=1:(N-1)
                NUM=N-I1;
                PA1A2=pa1a2(A1,A2,FA,LAM,I1,0);
                COM=NUM*2*PA1A2*Residue(I)*(LAMI^(1-G)*(1-LAMI^G)^2/(1-LAMI)^2)*LAMI^(I1*G);
                COM(isnan(COM))=0;
                val(j)=val(j)+COM;
            end
        end
    end
    val(j)=val(j)/(N*G);
    
    if abs(imag(val(j)))<1e-14
        val(j)=real(val(j));
    end
end
