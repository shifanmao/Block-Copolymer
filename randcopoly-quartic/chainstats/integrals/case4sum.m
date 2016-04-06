function valne=case4sum(N,Z3,Z3L)
    if N<2
        valne=zeros(2,2,2);
        return
    end
    valne1=twoSum(Z3,N);
    valne2=twoSum(Z3L,N);
    valne=zeros(2,2,2);
    valne(:,:,1)=ones(2,2)*valne1;
    valne(:,:,2)=ones(2,2)*valne2;
end