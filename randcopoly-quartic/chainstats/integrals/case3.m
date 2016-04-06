function S4=case3(N,NM,R1,R12,R3,Z12,Z12L,F,order)
% Case 3: J1==J2 <J3==J4
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE12=Z12;
    ZE12L=Z12L;

    valeq=case3Int(E1,E12,E3,NM);
    
    % on different monomers
    valne1=twoSum(ZE12,N);
    valne2=twoSum(ZE12L,N);
    
    valne=zeros(2,2,2);
    valne(:,1,:)=ones(2,2)*valne1;
    valne(:,2,:)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end
