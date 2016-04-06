function s2=s2wlc(N,NM,LAM,FA,k)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

% calculate the eigenvalues
if sum(power(k,2)) < 1e-8
    k=1e-5;
end

ORDEig=20;  % maximum number of eigenvalues
ResLayer=500;  % number of residual layers
R=Eigenvalues(k,ORDEig,1);
NR=ORDEig;

% get the residues for all roots of each k(j)
Residue=Residues(k,R(1:NR),NR,1,ResLayer,1);

FB=1-FA;
F=[FA,FB];
DF=[1,-1];
s2=zeros(2,2);
for I=1:NR
    Z0=exp(R(I)*NM);
    Z0L=Z0*LAM;        
% Case 1 :: J1==J2
    % on same monomer (S integrals)
    valeq = R(I).^(-2).*expl(2,R(I)*NM);

    % on different monomers (J sums)
    valne = N;

    s2(1,1)=s2(1,1)+2*Residue(I)*F(1)*valeq*valne;
    s2(2,2)=s2(2,2)+2*Residue(I)*F(2)*valeq*valne;

% Case 2 :: J1~=J2        
    if (NM*R(I)) < -10
        valeq=(Z0-1)^2.*R(I)^(-2);
        valne1=(-1).*((-1)+Z0).^(-2).*Z0.*(1+N.*((-1)+Z0)+(-1).*Z0.^N);
        valne2=(-1).*((-1)+Z0L).^(-2).*LAM.*(1+N.*((-1)+Z0L)+(-1).*Z0L.^N);
    else 
        % on same monomer (S integrals)
        valeq = R(I).^(-2).*(expl(2,R(I)*NM)+expl(2,-R(I)*NM));

        % on different monomers (J sums)
        if  Z0>0.7
            % Use form for Z0 near 1, note Z0 <= 1 always
            valne1 = (-Z0).*(Z0-1).^(-2).*(-expl(2,N.*R(I)*NM)+expl(2,R(I)*NM).*N);            
        else
            % Use form for Z0 closer to 0 to avoid log(Z0)
            valne1 = (-1).*((-1)+Z0).^(-2).*Z0.*(1+N.*((-1)+Z0)+(-1).*Z0.^N);
        end

        if Z0L>0.7
            % Use form for Z0L near 1, note Z0L <= 1 always
            valne2 = (-Z0L).*(Z0L-1).^(-2).*(-expl(2,N.*log(Z0L))+expl(2,log(Z0L)).*N);            
        else
            % Use form for Z0L closer to 0 or negitive to avoid log(Z0)
            valne2 = (-1).*((-1)+Z0L).^(-2).*Z0L.*(1+N.*((-1)+Z0L)+(-1).*Z0L.^N);
        end
    end
    if isnan(Residue(I))
        error('Residue is nan')
    end
    if isnan(R(I))
        error('R is nan')
    end
    if isnan(valeq)
        error('valeq is nan')
    end
    if isnan(valne1)
        sprintf('Value of Z0 %f',Z0)
        error('valne1 is nan')
    end
    if isnan(valne2)
        error('valne2 is nan')
    end

    % J1<J2
    s2(1,1)=s2(1,1)+Residue(I)*(F(1)*F(1)*valeq*valne1+F(1)*DF(1)*DF(1)*(1-F(1))*valeq*valne2);
    s2(2,2)=s2(2,2)+Residue(I)*(F(2)*F(2)*valeq*valne1+F(2)*DF(2)*DF(2)*(1-F(2))*valeq*valne2);
    s2(1,2)=s2(1,2)+Residue(I)*(F(1)*F(2)*valeq*valne1+F(1)*DF(1)*DF(2)*(1-F(1))*valeq*valne2);

    %J2<J1
    s2(1,1)=s2(1,1)+Residue(I)*(F(1)*F(1)*valeq*valne1+F(1)*DF(1)*DF(1)*(1-F(1))*valeq*valne2);
    s2(2,2)=s2(2,2)+Residue(I)*(F(2)*F(2)*valeq*valne1+F(2)*DF(2)*DF(2)*(1-F(2))*valeq*valne2);
    s2(1,2)=s2(1,2)+Residue(I)*(F(2)*F(1)*valeq*valne1+F(2)*DF(2)*DF(1)*(1-F(2))*valeq*valne2); 

end
s2(2,1)=s2(1,2);

if(sum(sum(isnan(s2))))
    format longG
    sprintf('valeq=%f, valne1=%f, valne2=%f, Z0=%f, Z0L=%f, R(I)=%f',...
            valeq,valne1,valne2,Z0,Z0L,R(I))
    sprintf('s2:')
    disp(s2);
    error('s2 is nan');
end
    
end