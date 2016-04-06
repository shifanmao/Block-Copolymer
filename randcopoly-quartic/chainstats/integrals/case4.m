function S4=case4(N,NM,R1,R12,R3,Z3,Z3L,F,order)
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    % Case 4: J1==J2==J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;
    ZE3=Z3;
    ZE3L=Z3L;
   
    valeq=case4Int(E1,E12,E3,NM);
    % on different monomers
    valne1=twoSum(ZE3,N);
    valne2=twoSum(ZE3L,N);

    valne=zeros(2,2,2);
    valne(:,:,1)=ones(2,2)*valne1;
    valne(:,:,2)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end