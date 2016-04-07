function S4=s4wlc(N,NM,LAM,FA,Q1,Q2,Q3,Q4)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA
ORDEig=4;  % maximum number of eigenvalues
ORDL=4;
NumLayer=500;  % number of residual layers

% Reset Qs to column vectors if entered as rows
if isarow(Q1)==1
    Q1=transpose(Q1);
    Q2=transpose(Q2);
    Q3=transpose(Q3);
    Q4=transpose(Q4);
end

% Begin calculation of s4

% exclude zero vector
MIN=1e-6;
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3+Q4,2)))])
    error('Wavevectors must add up to zero from translational invariance')
end

% Evaluate the quantities for s4 calculation
FB=1-FA;
F=[FA,FB];

S4=zeros(2,2,2,2);
%S4split=zeros(8,2,2,2,2);
orders = perms(1:4);
for orderNum=1:24
    
    order=orders(orderNum,:);
    Q = [Q1,Q2,Q3,Q4];
    Qnew=[Q(:,order(1)),Q(:,order(2)),Q(:,order(3)),Q(:,order(4))];
    % Qnew is the reordered Q

    % Now take advantage of translatioanl invariance
    q1=-Qnew(:,1);
    q2=-(Qnew(:,1)+Qnew(:,2));
    q3=Qnew(:,4);

    % Find unit vectors
    Q1_n=q1/norm(q1,2);
    Q2_n=q2/norm(q2,2);
    Q3_n=q3/norm(q3,2);

    % Find Euler angles
    aligned=10^-13; % angle in radians betwene to vectors before I assume they are the same
    if norm(q2) < aligned
        % if q2=0 we may as well point it perpendicular to the other two
        cosB1=0;
        cosB2=0;
        alpha1=dot(Q1_n,Q3_n);
    elseif abs(Q2_n-Q1_n) < aligned
        cosB1=1;
        alpha1=0;
        cosB2=dot(Q2_n,Q3_n);
    elseif abs(Q2_n-Q3_n) < aligned
        cosB1=dot(Q2_n,Q1_n);
        alpha1=0;
        cosB2=1;
    else
        cosB1=dot(Q2_n,Q1_n);
        cosB2=dot(Q2_n,Q3_n);
        v1=cross(Q2_n,Q1_n)/norm(cross(Q2_n,Q1_n));
        v2=cross(Q2_n,Q3_n)/norm(cross(Q2_n,Q3_n));
        alpha1=acos(dot(v1,v2));
    end 
    % Calculate Wigner D matrices
    firstD=WignerD_lm0(ORDL,alpha1,cosB1);
    secondD=WignerD_lm0(ORDL,0,cosB2);
    % index: (lam+1,M+1)
    
    % Now calculate the eigenvalue
    R1mtrx=Eigenvalues(norm(q1),ORDEig,ORDL); % won't need all of output
    R12mtrx=Eigenvalues(norm(q2),ORDEig,ORDL);
    R4mtrx=Eigenvalues(norm(q3),ORDEig,ORDL); % won't need all of output
    % index: (l+1,M+1)
    
    GL1=Residues(norm(q1),R1mtrx,ORDEig,ORDL,NumLayer,1); 
    GLM12=Residues(norm(q2),R12mtrx,ORDEig,ORDL,NumLayer,ORDL);
    GL4=Residues(norm(q3),R4mtrx,ORDEig,ORDL,NumLayer,1);  
    % index: (lam1+1, lam2+1, mu+1, l+1)
    
    if mod(ORDEig,2)==1
        error('please just make my life easy and make ORDEig even!')
    end
    for M=0:(ORDL-1)
        if mod(ORDEig-M,2)==0
            Lmax=ORDEig-1;
        else
            Lmax=ORDEig-2;
        end
        
        for L1=0:(ORDEig-1)
            for L2=M:Lmax
                for L3=0:(ORDEig-1) 
                    
                    R1=  R1mtrx(L1+1,1);
                    R12=R12mtrx(L2+1,M+1);
                    R4=  R4mtrx(L3+1,1);
                    
                    Z1= exp(R1*NM);
                    Z12=exp(R12*NM);
                    Z4= exp(R4*NM);
                    Z1L=  Z1*LAM;
                    Z12L=Z12*LAM;
                    Z4L=  Z4*LAM;
                    
                    if min(real([R1,R12,R4]*NM))<-700
                        error('overflow warning, need to cut off small epsilon')
                    end

                    valsne=zeros(2,2,2,8);
                    valseq=zeros(8,1);
                    
                    valseq(1)=case1Int(R1,R12,R4,NM);
                    valsne(:,:,:,1)=ones(2,2,2)*N;
                    
                    
                    if N>1; 
                    valseq(2)=case4Int(R4,R12,R1,NM);
                    valsne(:,:,:,2)=case2sum(N,Z1,Z1L);
                    
                    valseq(3)=case3Int(R1,R12,R4,NM);
                    valsne(:,:,:,3)=case3sum(N,Z12,Z12L);
                    
                    valseq(4)=case4Int(R1,R12,R4,NM); % Case 4: J1==J2==J3 <J4
                    valsne(:,:,:,4)=case4sum(N,Z4,Z4L);
                    end
                    if N>2
                    valseq(5)=case6Int(R4,R12,R1,NM); % Case 5: J1 <J2 <J3==J4
                    valsne(:,:,:,5)=case5sum(N,Z1,Z1L,Z12,Z12L);
                    
                    valseq(6)=case6Int(R1,R12,R4,NM); % Case 6: J1==J2 <J3 <J4
                    valsne(:,:,:,6)=case6sum(N,Z12,Z12L,Z4,Z4L);
                    
                    valseq(7)=case7Int(R1,R12,R4,NM);% Case 7: J1 <J2==J3 <J4
                    valsne(:,:,:,7)=case7sum(N,Z1,Z1L,Z4,Z4L);
                    end
                    if N>3
                    valseq(8)=case8Int(R1,R12,R4,NM);
                    valsne(:,:,:,8)=case8sum(N,Z1,Z1L,Z12,Z12L,Z4,Z4L);
                    end
                    S4term=BinomialSum(valseq,valsne,order,F);
                    
                    if any(any(any(any(any(isnan(S4term))))))
                        disp('valseq')
                        disp(valseq)
                        disp('S4term')
                        disp(S4term)
                        error('S4 term is NaN')
                    end
                    RotTerm=0;
                    if M==0 % If statment needed to include negitive \mu values
                        for lam2=M:(ORDL-1)
                            for lam3=M:(ORDL-1)
                                RotTerm=RotTerm+firstD(lam2+1,1)...
                                               *conj(secondD(lam3+1,1))...
                                               *GL1(1,lam2+1,1,L1+1)...
                                               *GLM12(lam2+1,lam3+1,1,L2+1)...
                                               *GL4(lam3+1,1,1,L3+1);
                                if isnan(RotTerm)
%                                     Q1_n
%                                     Q2_n
%                                     Q3_n
%                                     alpha1
%                                     cosB1
%                                     cosB2
%                                     firstD
%                                     GL1
%                                     GLM12
%                                     GL4
                                    error('NaN encountered')
                                end
                            end
                        end
                    else
                        for lam2=M:(ORDL-1)
                            for lam3=M:(ORDL-1)
                                RotTerm=RotTerm+2*real(firstD(lam2+1,M+1)...
                                               *conj(secondD(lam3+1,M+1)) )...
                                               *GL1(1,lam2+1,1,L1+1)...
                                               *GLM12(lam2+1,lam3+1,M+1,L2+1)...
                                               *GL4(lam3+1,1,1,L3+1);
                                           
                                if isnan(RotTerm)
                                    error('NaN encountered')
                                end
                            end
                        end
                    end
                    
                    S4=S4+sum(S4term,5)*RotTerm;
                end
            end
        end
    end
   
end

end

function out=WignerD_lm0(ORDL,alpha,cosB)
% Wigner D Matrix when second L index is zero
% Returns a ORDL x ORDL matrix with indicies (lam+1,M+1)
out=zeros(ORDL,ORDL)*NaN;

for L=0:(ORDL-1)
    M=(0:L);
    out(L+1,1:(L+1))=legendre(L,cosB)'.*...
                     sqrt(factorial(L-M)./factorial(L+M))...
                     .*exp(1i*M*alpha);
end
end