function valne=case6sum(N,Z12,Z12L,Z3,Z3L)
    if N<3
        valne=zeros(2,2,2);
        return
    end
    % Case 6: J1==J2 <J3 <J4          
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    
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
end