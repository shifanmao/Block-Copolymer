function S4=case5(N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,F,order)
    if N<3
        S4=zeros(2,2,2,2);
        return
    end
    % Case 5: J1 <J2 <J3==J4 
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];


    valeq=case6Int(E3,E12,E1,NM);
  
    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT12=ZE12(I2);
            valne(I1,I2,1)=tripleSum(ZT1,ZT12,N);
            valne(I1,I2,2)=valne(I1,I2,1);
        end
    end
    
    S4=BinomialSum(valeq,valne,order,F);
    
end