clc


NM=0.1;
vals=[0,-20,-1.1*10^-4/NM,-0.9*10^-4/NM,-0.5+0.4*10^-4/NM,-0.5,-0.5-0.4*10^-4/NM,-0.6];


R1=vals(1);
R12=vals(1);
R4=vals(3);



N=5;

a=exp(R1*NM);
b=exp(R12*NM);
c=exp(R4*NM);






val=zeros(8,2);

val(1,1)=case1Int(R1,R12,R4,NM);
val(1,2)=N;

if N>1; 
val(2,1)=case4Int(R4,R12,R1,NM);
val(2,2)=twoSum(a,N);

val(3,1)=case3Int(R1,R12,R4,NM);
val(3,2)=twoSum(b,N);

val(4,1)=case4Int(R1,R12,R4,NM); % Case 4: J1==J2==J3 <J4
val(4,2)=twoSum(c,N);
end
if N>2
val(5,1)=case6Int(R4,R12,R1,NM); % Case 5: J1 <J2 <J3==J4
val(5,2)=tripleSum(a,b,N);

val(6,1)=case6Int(R1,R12,R4,NM); % Case 6: J1==J2 <J3 <J4
val(6,2)=tripleSum(b,c,N);

val(7,1)=case7Int(R1,R12,R4,NM);% Case 7: J1 <J2==J3 <J4
val(7,2)=tripleSum(a,c,N);
end
if N>3
val(8,1)=case8Int(R1,R12,R4,NM);
val(8,2)=case8JPart(a,b,c,N);
end
format longEng
disp('all')
disp(val)
disp('sums')
disp(sum(val))
disp('prod')
disp(val(:,1)'*val(:,2))