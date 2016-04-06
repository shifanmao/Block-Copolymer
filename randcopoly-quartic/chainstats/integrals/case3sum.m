function valne=case3sum(N,Z12,Z12L)
    if N<2
        valne=zeros(2,2,2);
        return
    end
    valne1=twoSum(Z12,N);
    valne2=twoSum(Z12L,N);
    valne=zeros(2,2,2);
    valne(:,1,:)=ones(2,2)*valne1;
    valne(:,2,:)=ones(2,2)*valne2;
end