function s2=s2gc(N,NM,LAM,FA,k)

% Calculate the Fourier transform of the Green function
% for the gaussian chain in d-dimension
%
% Andrew Spakowitz (4/14/15)
s2=zeros(2,2);
d=3;    % three dimensional
R =-k*k/(2*d);
Z0=exp(R*NM);
Z0L=Z0*LAM;

FB=1-FA;
F=[FA,FB];
DF=[1,-1];

% Case 1 :: J1==J2
    % on same monomer (S integrals)
    %valeq = Z0/R^2-1/R^2-NM/R;% comment QJM 6/22/15
    valeq = R.^(-2).*expl(2,R*NM); %QJM 6/22/15

    % on different monomers (J sums)
    valne = N;

    s2(1,1)=s2(1,1)+2*F(1)*valeq*valne;
    s2(2,2)=s2(2,2)+2*F(2)*valeq*valne;

% Case 2 :: J1~=J2        
    if (NM*R) < -10
        % canceled evaluation for limit for large negitive R(I)*NM
        % Don't need to worry about z near 1 numerical issue
        valeq=(Z0-1)^2.*R^(-2);
        valne1=(-1).*((-1)+Z0).^(-2).*Z0.*(1+N.*((-1)+Z0)+(-1).*Z0.^N);
        valne2=(-1).*((-1)+Z0L).^(-2).*LAM.*(1+N.*((-1)+Z0L)+(-1).*Z0L.^N);
    else 
        % on same monomer (S integrals)
        %valeq = R.^(-2).*((-1)+Z0).^2.*Z0.^(-1); % comment QJM 6/22/15
        valeq = R^(-2).*(expl(2,R*NM)+expl(2,-R*NM));

        % on different monomers (J sums)
        if  Z0>0.7
            % Use form for Z0 near 1, note Z0 <= 1 always
            valne1 = (-Z0).*(Z0-1).^(-2).*(-expl(2,N.*R*NM)+expl(2,R*NM).*N);            
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


    % J1<J2
    s2(1,1)=s2(1,1)+F(1)*F(1)*valeq*valne1+F(1)*DF(1)*DF(1)*(1-F(1))*valeq*valne2;
    s2(2,2)=s2(2,2)+F(2)*F(2)*valeq*valne1+F(2)*DF(2)*DF(2)*(1-F(2))*valeq*valne2;
    s2(1,2)=s2(1,2)+F(1)*F(2)*valeq*valne1+F(1)*DF(1)*DF(2)*(1-F(1))*valeq*valne2;
    
    %J2<J1
    s2(1,1)=s2(1,1)+F(1)*F(1)*valeq*valne1+F(1)*DF(1)*DF(1)*(1-F(1))*valeq*valne2;
    s2(2,2)=s2(2,2)+F(2)*F(2)*valeq*valne1+F(2)*DF(2)*DF(2)*(1-F(2))*valeq*valne2;
    s2(1,2)=s2(1,2)+F(2)*F(1)*valeq*valne1+F(2)*DF(2)*DF(1)*(1-F(2))*valeq*valne2;
    
s2(2,1)=s2(1,2);

end
