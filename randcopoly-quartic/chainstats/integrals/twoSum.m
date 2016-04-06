function out=twoSum(a,N)
   if N<2
       out=0;
       return
   elseif abs(1-a)>0.2
       out=a*(N-N*a+a^N-1)/((1-a)^2);
   elseif abs(1-a)<10^-13
       out=N*(N+1)/2 -N;
   else
       out=a*(-N*expl(2,log(a))+expl(2,N*log(a)))/((a-1)^2);
   end
end