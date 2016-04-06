function S4=case6(N,NM,R1,R12,R3,Z12,Z12L,Z3,Z3L,F,order)
    if N<3
        S4=zeros(2,2,2,2);
        return
    end
    % Case 6: J1==J2 <J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
   
    % on same monomer\
    valeq=case6Int(E1,E12,E3,NM);
    
    % on different monomers
    valne=zeros(2,2,2);
    for I2=1:2
        for I3=1:2
            ZT12=ZE12(I2);
            ZT3=ZE3(I3);
            valne(1,I2,I3)=tripleSum(ZT12,ZT3,N);
            valne(2,I2,I3)=valne(1,I2,I3);
        end
    end
%     valeq=finite(valeq);
%     valne=finite(valne);

    S4=BinomialSum(valeq,valne,order,F);
end