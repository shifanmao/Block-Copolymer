function valne=case2sum(N,Z1,Z1L)
    if N<2
        valne=zeros(2,2,2);
        return
    end
    valne1=twoSum(Z1,N);
    valne2=twoSum(Z1L,N);  
    valne=zeros(2,2,2);
    valne(1,:,:)=ones(2,2)*valne1;
    valne(2,:,:)=ones(2,2)*valne2;
end