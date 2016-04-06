function S4=s4gc(N,NM,LAM,FA,Q1,Q2,Q3,Q4)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

d=3;    % three dimensional

% Begin calculation of s4
MIN=1e-6;
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3+Q4,2)))])
    error('Wavevectors must add up to zero from translational invariance')
end

% Evaluate the quantities for s4 calculation
FB=1-FA;
F=[FA,FB];

S4=zeros(2,2,2,2);
orders = perms(1:4);
for orderNum=1:24
    order=orders(orderNum,:);
    Q = [Q1,Q2,Q3,Q4];
    Qnew=[Q(:,order(1)),Q(:,order(2)),Q(:,order(3)),Q(:,order(4))];
    % Qnew is the reordered Q

    % Now calculate the eigenvalues 
    R1=-dot(Qnew(:,1),Qnew(:,1))/(2*d);
    R12=-dot(Qnew(:,1)+Qnew(:,2),Qnew(:,1)+Qnew(:,2))/(2*d);
    R4=-dot(Qnew(:,4),Qnew(:,4))/(2*d);

    Z1=exp(R1*NM);
    Z12=exp(R12*NM);
    Z4=exp(R4*NM);
    Z1L=Z1*LAM;
    Z12L=Z12*LAM;
    Z4L=Z4*LAM;

    S4=S4+case1(N,NM,R1,R12,R4,F);
    S4=S4+case2(N,NM,R1,R12,R4,Z1,Z1L,                F,order);
    S4=S4+case3(N,NM,R1,R12,R4,       Z12,Z12L,       F,order);
    S4=S4+case4(N,NM,R1,R12,R4,                Z4,Z4L,F,order);
    S4=S4+case5(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,       F,order);
    S4=S4+case6(N,NM,R1,R12,R4,       Z12,Z12L,Z4,Z4L,F,order);
    S4=S4+case7(N,NM,R1,R12,R4,Z1,Z1L,         Z4,Z4L,F,order);
    S4=S4+case8(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,order);
end
    
end
