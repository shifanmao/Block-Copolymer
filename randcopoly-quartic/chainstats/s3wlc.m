function S3=s3wlc(N,NM,LAM,FA,Q1,Q2,Q3)
S3=zeros(2,2,2);
MIN=5e-4;

ORDEig=20;  % maximum number of eigenvalues
ORDL=20;
NumLayer=500;  % number of residual layers

% Begin calculation of s3
if sum(power(Q1+Q2+Q3,2)) > MIN
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3,2)))])
    error('Wavevectors must add up to zero from translational invariance')
else
    if sum(power(Q3,2)) < MIN
        Q3=Q1*1e-5;
    end
    
    % Evaluate the quantities for s3 calculation
    FB=1-FA;
    F=[FA,FB];
    DF=[1,-1];
    
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));

    EQ1=Q1/Q1MAG;
    EQ2=Q2/Q2MAG;
    EQ3=Q3/Q3MAG;
    
    RHO12=sum(EQ1.*EQ2);
    RHO23=sum(EQ2.*EQ3);
    RHO13=sum(EQ1.*EQ3);
    
    PL12=legendrep(-RHO12,ORDL); % legendre Polynomial
    PL23=legendrep(-RHO23,ORDL); % P(n,x)
    PL13=legendrep(-RHO13,ORDL);

    % calculate the eigenvalues
    R1=Eigenvalues(Q1MAG,ORDEig,1);
    R2=Eigenvalues(Q2MAG,ORDEig,1);
    R3=Eigenvalues(Q3MAG,ORDEig,1);
    
    % get the residues for all roots of each k(j)
    GL1=Residues(Q1MAG,R1,ORDEig,ORDL,NumLayer,1);
    GL2=Residues(Q2MAG,R2,ORDEig,ORDL,NumLayer,1);
    GL3=Residues(Q3MAG,R3,ORDEig,ORDL,NumLayer,1);
    
    % but only need one -> m=0
    GL1=GL1(:,1,1,:);
    GL2=GL2(:,1,1,:);
    GL3=GL3(:,1,1,:);
    
    Z1=exp(R1*NM);
    Z2=exp(R2*NM);
    Z3=exp(R3*NM);
    
    Z1L=exp(R1*NM)*LAM;
    Z2L=exp(R2*NM)*LAM;
    Z3L=exp(R3*NM)*LAM;
    
    for N1=1:ORDEig
        for N2=1:ORDEig
            % Case 1: J1=J2=J3
            S3=case1(S3,N,NM,N1,N2,R1,R2,PL12,GL1,GL2,F,MIN);
            S3=case1(S3,N,NM,N1,N2,R1,R3,PL13,GL1,GL3,F,MIN);
            S3=case1(S3,N,NM,N1,N2,R2,R3,PL23,GL2,GL3,F,MIN);

            if any(isnan(S3))
                error('Case 1:: s3 is nan')
            end
            
            % Case 2: J1<J2=J3
            % J1<J2=J3
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(2,1,1)=1;
            SDEL(1,2,2)=1;
            S3=case2(S3,N,NM,N1,N2,R1,R2,Z1,Z1L,PL12,GL1,GL2,F,DF,SDEL,1,2,MIN);
            S3=case2(S3,N,NM,N1,N2,R1,R3,Z1,Z1L,PL13,GL1,GL3,F,DF,SDEL,1,2,MIN);

            % J2<J1=J3
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(2,1,2)=1;
            SDEL(1,2,1)=1;
            S3=case2(S3,N,NM,N1,N2,R2,R3,Z2,Z2L,PL23,GL2,GL3,F,DF,SDEL,1,2,MIN);
            S3=case2(S3,N,NM,N1,N2,R2,R1,Z2,Z2L,PL12,GL2,GL1,F,DF,SDEL,1,2,MIN);

            % J3<J1=J2
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(1,1,2)=1;
            SDEL(2,2,1)=1;
            S3=case2(S3,N,NM,N1,N2,R3,R1,Z3,Z3L,PL13,GL3,GL1,F,DF,SDEL,3,1,MIN);
            S3=case2(S3,N,NM,N1,N2,R3,R2,Z3,Z3L,PL23,GL3,GL2,F,DF,SDEL,3,1,MIN);
            
            if any(isnan(S3))
                error('Case 2:: s3 is nan')
            end

            % Case 3: J1=J2<J3
            % J1=J2<J3
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(1,1,2)=1;
            SDEL(2,2,1)=1;
            S3=case3(S3,N,NM,N1,N2,R1,R3,Z3,Z3L,PL13,GL1,GL3,F,DF,SDEL,2,3,MIN);
            S3=case3(S3,N,NM,N1,N2,R2,R3,Z3,Z3L,PL23,GL2,GL3,F,DF,SDEL,2,3,MIN);

            % J1=J3<J2
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(1,2,1)=1;
            SDEL(2,1,2)=1;
            S3=case3(S3,N,NM,N1,N2,R1,R2,Z2,Z2L,PL12,GL1,GL2,F,DF,SDEL,3,2,MIN);
            S3=case3(S3,N,NM,N1,N2,R3,R2,Z2,Z2L,PL23,GL3,GL2,F,DF,SDEL,3,2,MIN);

            % J2=J3<J1
            SDEL=zeros(2,2,2);
            SDEL(1,1,1)=1;
            SDEL(2,2,2)=1;
            SDEL(1,2,2)=1;
            SDEL(2,1,1)=1;
            S3=case3(S3,N,NM,N1,N2,R2,R1,Z1,Z1L,PL12,GL1,GL2,F,DF,SDEL,3,1,MIN);
            S3=case3(S3,N,NM,N1,N2,R3,R1,Z1,Z1L,PL23,GL2,GL3,F,DF,SDEL,3,1,MIN);
            
            if any(isnan(S3))
                error('Case 3:: s3 is nan')
            end

            % Case 4: J1<J2<J3
            SDEL=ones(2,2,2);
            % J1<J2<J3
            S3=case4(S3,N,NM,N1,N2,R1,R3,Z1,Z3,Z1L,Z3L,PL13,GL1,GL3,F,DF,SDEL,1,2,3,MIN);

            %J3<J2<J1
            S3=case4(S3,N,NM,N1,N2,R3,R1,Z3,Z1,Z3L,Z1L,PL13,GL3,GL1,F,DF,SDEL,3,2,1,MIN);

            %J1<J3<J2
            S3=case4(S3,N,NM,N1,N2,R1,R2,Z1,Z2,Z1L,Z2L,PL12,GL1,GL2,F,DF,SDEL,1,3,2,MIN);

            %J2<J3<J1
            S3=case4(S3,N,NM,N1,N2,R2,R1,Z2,Z1,Z2L,Z1L,PL12,GL2,GL1,F,DF,SDEL,2,3,1,MIN);

            %J2<J1<J3
            S3=case4(S3,N,NM,N1,N2,R2,R3,Z2,Z3,Z2L,Z3L,PL23,GL2,GL3,F,DF,SDEL,2,1,3,MIN);

            %J3<J1<J2
            S3=case4(S3,N,NM,N1,N2,R3,R2,Z3,Z2,Z3L,Z2L,PL23,GL3,GL2,F,DF,SDEL,3,1,2,MIN);
            
            if any(isnan(S3))
                error('Case 4:: s3 is nan')
            end
            
        end
    end
    S3=real(S3);
end
end

function S3=case1(S3,N,NM,N1,N2,E1,E2,PL12,GL1,GL2,F,MIN)
% Case 1: J1=J2=J3
    R1=E1(N1);
    R2=E2(N2);
    
    % on same monomer
    if abs(R1)<2*MIN && abs(R2)<2*MIN  % notice the factor of 2
        valeq=NM^3*(NM*R1+NM*R2+4)/24;
    elseif abs(R1-R2)<MIN
        valeq=(-2*expl(3,NM*R1)+NM*R1*expl(2,NM*R1))/(R1^3);
    elseif abs(R1)<MIN
        valeq=expl(3,NM*R2)/(R2^3);
    elseif abs(R2)<MIN
        valeq=expl(3,NM*R1)/(R1^3);
    else
        valeq=(R1^2*expl(2,NM*R2)-R2^2*expl(2,NM*R1))/(R1^2*R2^2*(R2-R1));
    end
    
    S3(1,1,1)=S3(1,1,1)+2*F(1)*N*valeq*sum(PL12.*GL1(:,N1).*GL2(:,N2));
    S3(2,2,2)=S3(2,2,2)+2*F(2)*N*valeq*sum(PL12.*GL1(:,N1).*GL2(:,N2));
end

function S3=case2(S3,N,NM,N1,N2,E1,E2,ZT1,ZT1L,PL12,GL1,GL2,F,DF,SDEL,A1,A2,MIN)
% Case 2: J1<J2=J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z1L=ZT1L(N1);
    
    % on same monomer
    if abs(R1)<2*MIN && abs(R2)<2*MIN  % notice the factor of 2
        valeq=NM^3*(2*NM*R2-NM*R1+6)/12;
    elseif abs(R1-R2)<MIN
        valeq=(NM*R1*expl(2,NM*R1)-2*coshl(3,NM*R1))/(R1^3);
    elseif abs(R1)<MIN
        valeq=NM*expl(2,NM*R2)/(R2^2);
    elseif abs(R2)<MIN
        valeq=(2*coshl(3,NM*R1)+NM*R1*expl(2,-NM*R1))/(R1^3);
    else
        valeq=(2*R2*coshl(2,NM*R1)+R1*expl(1,-NM*R1)*expl(1,NM*R2))/(R1^2*R2*(R1-R2));
    end
    
    jMIN=10^-7;
    ZE1=[Z1,Z1L];
    valne=[NaN,NaN];
    for I1=1:2
        ZT1=ZE1(I1);
        if abs(ZT1-1)<jMIN
            valne(I1)=N*(N-1)*((N+1)*ZT1-(N-2))/6;
        else
            valne(I1)=ZT1*powl(2,ZT1,N)/((ZT1-1)^2);
        end
    end
    valeq(isinf(valeq))=0;
    valeq(isnan(valeq))=0;
    valne(isinf(valne))=0;

    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                          valeq*(FV(A1)*FV(A2)*valne(1)+...
                             FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*valne(2))*...
                          sum(PL12.*GL1(:,N1).*GL2(:,N2));
            end
        end
    end
end

function S3=case3(S3,N,NM,N1,N2,E1,E2,ZT2,ZT2L,PL12,GL1,GL2,F,DF,SDEL,A2,A3,MIN)
% Case 3: J1=J2<J3
    R1=E1(N1);
    R2=E2(N2);
    Z2=ZT2(N2);
    Z2L=ZT2L(N2);
    
    % on same monomer
    if abs(R1)<2*MIN && abs(R2)<2*MIN  % notice the factor of 2
        valeq=NM^3*(2*NM*R1-NM*R2+6)/12;
    elseif abs(R1-R2)<MIN
        valeq=(NM*R1*expl(2,NM*R1)-2*coshl(3,NM*R1))/(R1^3);
    elseif abs(R1)<MIN
        valeq=-expl(1,-NM*R2)*expl(2,NM*R2)/(R2^3);
    elseif abs(R2)<MIN
        valeq=NM*expl(2,NM*R1)/(R1^2);
    else
        valeq=(2*R1*coshl(2,NM*R2)+R2*expl(1,NM*R1)*expl(1,-NM*R2))/(R1*R2^2*(R2-R1));
    end
    
    % on different monomers
    jMIN=10^-7;
    ZE2=[Z2,Z2L];
    valne=[NaN,NaN];
    for I2=1:2
        ZT2=ZE2(I2);
        if abs(ZT2-1)<jMIN
            valne(I2)=N*(N-1)*((N+1)*ZT2-(N-2))/6;
        else
            valne(I2)=ZT2*powl(2,ZT2,N)/((ZT2-1)^2);
        end
    end
    valeq(isinf(valeq))=0;
    valeq(isnan(valeq))=0;
    valne(isinf(valne))=0;
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                          valeq*(FV(A2)*FV(A3)*valne(1)+...
                             FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(2))*...
                          sum(PL12.*GL1(:,N1).*GL2(:,N2));
            end
        end
    end
end

function S3=case4(S3,N,NM,N1,N2,E1,E2,ZT1,ZT2,ZT1L,ZT2L,PL12,GL1,GL2,F,DF,SDEL,A1,A2,A3,MIN)
% Case 4: J1<J2<J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z2=ZT2(N2);
    Z1L=ZT1L(N1);
    Z2L=ZT2L(N2);
    
    ZE1=[Z1,Z1L];
	ZE2=[Z2,Z2L];
    
    % on same monomer
    if abs(R1)<MIN && abs(R2)<MIN  % notice the factor of 2
        valeq=NM^3;
    elseif abs(R1)<MIN
        valeq=2*NM*coshl(2,NM*R2)/(R2^2);
    elseif abs(R2)<MIN
        valeq=2*NM*coshl(2,NM*R1)/(R1^2);
    elseif abs(R1-R2)<MIN
        valeq=2*NM*coshl(2,NM*R1)/(R1^2);
    else
        valeq=expl(1,NM*R1)*expl(1,-NM*R2)*expl(1,NM*(R2-R1))/(R1*R2*(R1-R2));
    end
    if (valeq>1e20)
        valeq=0;
    end

    % on different monomers
    valne=zeros(2,2)*NaN;
    %valne(1,1) = no lambda term
    %valne(1,2) = lambda^(J3-J2) term
    jMIN=10^-7;
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT2=ZE2(I2);

            if abs(ZT1-1)<2*jMIN && abs(ZT2-1)<2*jMIN
                valne(I1,I2)=N*(N-1)*(N-2)*((N+1)*(ZT1+ZT2)-2*(N-1))/24;
            elseif abs(ZT1-1)<jMIN
                valne(I1,I2)=ZT2*powl(3,ZT2,N)/((ZT2-1)^3);
            elseif abs(ZT2-1)<jMIN
                valne(I1,I2)=ZT1*powl(3,ZT1,N)/((ZT1-1)^3);
            elseif abs(ZT1-ZT2)<jMIN
                valne(I1,I2)=-(N*(N-1)*ZT1*(ZT1-1)^3+2*ZT1^2*powl(3,ZT1,N)-N*ZT1*(ZT1-1)*powl(2,ZT1,N))/((ZT1-1)^3);
            else
                valne(I1,I2)=-ZT1*ZT2*((ZT1-1)^2*powl(3,ZT2,N)-(ZT2-1)^2*powl(3,ZT1,N))/((ZT1-ZT2)*(ZT1-1)^2*(ZT2-1)^2);
            end
        end
    end
    valeq(isinf(valeq))=0;
    valne(isinf(valne))=0;
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                   valeq*(FV(A1)*FV(A2)*FV(A3)*valne(1,1)+...
                      FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(1,2)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*valne(2,1)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(2,2))*...
                   sum(PL12.*GL1(:,N1).*GL2(:,N2));
            end
        end
    end
end