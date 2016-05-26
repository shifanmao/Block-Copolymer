function Res=NewCalRes(L_inf,R,k)

% calculate the residual of all eigenvalues

ResThreshold=1e-12;     % threshold to go from small k asymptot to matrix method
NR=length(R);           % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative
d=3;                % dimensionality of the problem
Res=zeros(NR,1);    % residual vector, each row corresponds to each eigenvalue
NumLayers=L_inf;    % number of layers in the continued franction

% find the residual
for i=1:NR

    % use asymptotic residual to figure out whether calculate the residual
    % using continued fraction or stay with the asymptotic form for small k
    % limit
    Res(i)=SmallAsympRes(k,i);
    if abs(Res(i))>ResThreshold
        if k<=1

            % residual using continued fraction for small k
            p=R(i);
            W=p+(L_inf+d-2)*L_inf;
            Wprime=1;
            for L=NumLayers:-1:1
                AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
                Wprime=1-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm+AL^2/W;
            end
            Res(i)=1/Wprime;
        else

            % residual using continued fraction for large k
            p=R(i);
            W=(p+(L_inf+d-2)*L_inf)/k;
            Wprime=1/k;
            for L=NumLayers:-1:1
                AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
                Wprime=1/k-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm/k+AL^2/W;
            end
            Res(i)=1/(k*Wprime);
        end
    end
end

