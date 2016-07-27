function S3=s3wlc(NM,Q1,Q2,Q3,ORDEig,ORDL,ResLayer)
% Calculate the Fourier transform of the three-point Green function
% for the wormlike chain in d-dimension

S3=0;
NR=ORDEig;
MIN=1e-5;

% Begin calculation of s3
if sum(power(Q1+Q2+Q3,2)) > MIN
    
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3+Q4,2)))])
    error('Wavevectors must add up to zero from translational invariance')
    
elseif sum(power(Q3,2)) < MIN
    
    % reduces to two-point correlation function
    Q1MAG=sqrt(sum(power(Q1,2)));
    S3 = s2wlc(NM,Q1MAG,ORDEig,ResLayer)*NM;
    
else

    % Evaluate the quantities for s3 calculation
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));

    % calculate the eigenvalues
    R1=Eigenvalues(Q1MAG,NR,1);
    R2=Eigenvalues(Q2MAG,NR,1);
    R3=Eigenvalues(Q3MAG,NR,1);

    % get the residues for all roots of each k(j)
    GL1=Residues(Q1MAG,R1,ORDEig,ORDL,ResLayer,1);
    GL2=Residues(Q2MAG,R2,ORDEig,ORDL,ResLayer,1);
    GL3=Residues(Q3MAG,R3,ORDEig,ORDL,ResLayer,1);

    % but only need one -> m=0
    GL1=GL1(:,1,1,:);
    GL2=GL2(:,1,1,:);
    GL3=GL3(:,1,1,:);

    % unit vectors in direction of wavevectors
    EQ1=Q1/Q1MAG;
    EQ2=Q2/Q2MAG;
    EQ3=Q3/Q3MAG;

    % angles between wavevectors
    RHO12=sum(EQ1.*EQ2);
    RHO23=sum(EQ2.*EQ3);
    RHO13=sum(EQ1.*EQ3);

    % legendre polynomial representations of angles
    PL12=legendrep(-RHO12,ORDL);
    PL23=legendrep(-RHO23,ORDL);
    PL13=legendrep(-RHO13,ORDL);

    R1=NM*R1;
    R2=NM*R2;
    R3=NM*R3;
    for N1=1:NR
        for N2=1:NR
            % Case 1: A1=A2=A3 (SAAA)
            S3=S3_case1(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12);
            S3=S3_case1(S3,NM,N1,N2,R1,R3,GL1,GL3,PL13);
            S3=S3_case1(S3,NM,N1,N2,R2,R3,GL2,GL3,PL23);
            S3=S3_case1(S3,NM,N1,N2,R2,R1,GL2,GL1,PL12);
            S3=S3_case1(S3,NM,N1,N2,R3,R1,GL3,GL1,PL13);
            S3=S3_case1(S3,NM,N1,N2,R3,R2,GL3,GL2,PL23);
        end
    end

end
end

function S3=S3_case1(S3,NM,N1,N2,R1,R2,GL1,GL2,PL12)
% Case 1: A1=A2=A3
    S3=S3+case1_int(1,R1(N1),R2(N2))*NM^3*sum(PL12.*GL1(:,N1).*GL2(:,N2));
end

function valeq=case1_int(FA,R1,R2)
    MIN=(10^-2)/FA;
    if abs(R1-R2)<MIN
        if abs(R1)<MIN
            valeq=(1/6).*FA.^3;
        else
            valeq=R1.^(-3).*(2+FA.*R1+exp(1).^(FA.*R1).*((-2)+FA.*R1));
        end
    else
        if abs(R1)<MIN
            valeq=(-1/2).*R2.^(-3).*(2+(-2).*exp(1).^(FA.*R2)+2.*FA.*R2+FA.^2.*R2.^2);
        elseif abs(R2)<MIN
            valeq=(-1/2).*R1.^(-3).*(2+(-2).*exp(1).^(FA.*R1)+2.*FA.*R1+FA.^2.*R1.^2);
        else
            valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-2).*(((-1)+exp(1).^(FA.*R1)) ...
              .*R2.^2+(-1).*FA.*R1.*R2.^2+R1.^2.*(1+(-1).*exp(1).^(FA.*R2)+FA.*R2));
        end
    end
end