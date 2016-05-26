function Eig=MatRoots(k,N)

% find roots of denominator (eigenvalues) by solving eigenvalue problem

if k>8000

    % use large k asmyptotic expansion for large k
    d=3;                % dimensionality of the problem
    m=(d-3)/2;
    NumPoles=N;
    lambda=k;
    alpha=1/sqrt(8*lambda);
    I=complex(0,1);
    r=0;
   for np=1:NumPoles
        if r==NumPoles break; end
        r=r+1;
        n=2*(r-1)+m+1;
        Eig(r)=-(m*(m+1)*0-lambda*I+Epsilon(r-1,d,alpha));
        if imag(Eig(r))~=0
            r=r+1;
            Eig(r)=conj(Eig(r-1));
       end
    end
else

    % use matrix method for intermediate and small k regime
    n=4*N;
    E=zeros(n,n);
    for m=1:n
        if k<=1
            a=complex(0,-m*k/sqrt(4*m^2-1));
            if m>1 b=complex(0,-(m-1)*k/sqrt(4*(m-1)^2-1)); end
            if m==1
                E(m,1:2)=[m*(m-1),a];
            elseif m==n
                E(m,n-1:n)=[b,m*(m-1)];
            else
                E(m,m-1:m+1)=[b,m*(m-1),a];
            end
        else
            a=complex(0,-m/sqrt(4*m^2-1));
            if m>1 b=complex(0,-(m-1)/sqrt(4*(m-1)^2-1)); end
            if m==1
                E(m,1:2)=[m*(m-1)/k,a];
            elseif m==n
                E(m,n-1:n)=[b,m*(m-1)/k];
            else
                E(m,m-1:m+1)=[b,m*(m-1)/k,a];
            end
        end
    end
    TempMat=eig(E);
    [junk,index]=sort(real(TempMat));
    TempMat=TempMat(index);
    if k<=1
        Eig=-TempMat(1:N);
    else
        Eig=-TempMat(1:N)*k;
    end
end
