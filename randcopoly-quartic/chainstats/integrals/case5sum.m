function valne=case5sum(N,Z1,Z1L,Z12,Z12L)
    if N<3
        valne=zeros(2,2,2);
        return
    end
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];    
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT12=ZE12(I2);
            valne(I1,I2,1)=tripleSum(ZT1,ZT12,N);
            valne(I1,I2,2)=valne(I1,I2,1);
        end
    end
    
end