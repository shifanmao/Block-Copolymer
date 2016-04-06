function valne=case8sum(N,Z1,Z1L,Z12,Z12L,Z3,Z3L)
    if N<4
        valne=zeros(2,2,2);
        return
    end
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];

    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            for I3=1:2
                ZT1=ZE1(I1);
                ZT12=ZE12(I2);
                ZT3=ZE3(I3);
                valne(I1,I2,I3)=case8JPart(ZT1,ZT12,ZT3,N);
            end
        end
    end
end