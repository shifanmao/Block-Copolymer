function valne=case7sum(N,Z1,Z1L,Z3,Z3L)
    if N<3
        valne=zeros(2,2,2);
        return
    end
    % Case 7: J1 <J2==J3 <J4           
    ZE1=[Z1,Z1L];
    ZE3=[Z3,Z3L];
    
    valne=zeros(2,2,2);
    for I1=1:2
        for I3=1:2
            ZT1=ZE1(I1);
            ZT3=ZE3(I3);
            valne(I1,1,I3)=tripleSum(ZT1,ZT3,N);
            valne(I1,2,I3)=valne(I1,1,I3);
        end
    end
end