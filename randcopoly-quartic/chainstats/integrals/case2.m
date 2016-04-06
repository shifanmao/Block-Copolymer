function S4=case2(N,NM,R1,R12,R3,Z1,Z1L,F,order)
% Case 2: J1 <J2==J3==J4
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=Z1;
    ZE1L=Z1L;

    valeq=case4Int(E3,E12,E1,NM);
    % on different monomers
    valne1=twoSum(ZE1,N);
    valne2=twoSum(ZE1L,N);  
    
    valne=zeros(2,2,2);
    valne(1,:,:)=ones(2,2)*valne1;
    valne(2,:,:)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end
