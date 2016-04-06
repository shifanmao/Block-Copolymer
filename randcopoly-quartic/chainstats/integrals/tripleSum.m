function out=tripleSum(a,b,N)
%{
brute=0;
for j3=3:N
    for j2=2:j3-1
        for j1=1:j2-1
            brute=brute+(b^(j3-j2))*(a^(j2-j1));
        end
    end
end
%}

if abs(a*b)<1e-40
        out=0;
        return 
end

if N<3
    out=0;
    return
end

ndigits=14.5;  
if isreal(a) && isreal(b)
    ordered=sort([b,a]);
    b=ordered(1); 
    a=ordered(2);
else
    ab=[a,b];
    [junk,index]=min(abs(1-ab));
    a=ab(index);
    ab(index)=[];
    b=ab;
end

goodEnough=0.001;
% 1ab
oab=(-N*(N+1)*(N^2-3*N+2)/24)*(2-a-b)+...
    N*(N+1)*(4*N-16)/24+N;
oab_ac=max(N*(1-a),N*(1-b));
if oab_ac<goodEnough
    out=oab;
    return
end


% a~b~0
abz=-a*b*(2-N);
abz_ac=abs(a)+abs(b);
if abz_ac<goodEnough
    out=abz;
    return
end

% 1 !~ a !~ b~0
% o_a_bz=-a*b*(a^2*(N-1)-a*(N-2)-a^N)/((a-b)*(a-1)^2);
% o_a_bz_ac=

% 1 a b
M1=[...
    0   -1  2   -1;...
    1   0   -N  N-1;
    -2  N   0   2-N;...
    1   1-N N-2 0];
vL=[a^N,a^2,a,1];
vR=[b^N,b^2,b,1].';
prefac=-a*b/((a-b)*(1-a)^2*(1-b)^2);
o_a_b=vL*M1*vR*prefac;
ratio=abs( max(max(abs((vL.')*(vR.').*M1)))*prefac/o_a_b );
if isnan(ratio) || abs(a-b)<0.00001
    o_a_b_ac=100;
else
    o_a_b_ac=ratio*10^(-ndigits);
end
if o_a_b_ac<goodEnough
    out=o_a_b;
    return
end

% 1 ab
M1=[...
    0   0   0   2   -4;...
    -1  3   -2  -4  4;...
    2   -2  -4  2   0;...
    -1  -1  0   0   0;...
    0   0   0   2   0;...
    0   0   0   -4  4;...
    0   2   2   2   -4;...
    0   -2  4   0   0];
vR=[(a-b)*N^2, (a-b)*N, a-b,N,1].';
vL=(a*ones(1,8)).^([N+4,N+3,N+2,N+1,5,4,3,2]);
prefac=1/(2*a*(1-a)^4);
o_ab=vL*M1*vR*prefac;
ratio=abs(max(max(abs((vL.')*(vR.').*M1)))*prefac...
      /abs(o_ab));
if isnan(ratio)
        o_ab_ac=100;
else
    o_ab_ac=max(ratio*10^(-ndigits),N*abs((a-b)/a));
end
if o_ab_ac<goodEnough
    out=o_ab;
    return
end


% 1a b
M1=[...
    0   0   0   6   0   0   -6;...
    0   0   0   0   0   0   6;...
    -1  0   1   0   3   -3  0;...
    3   -3  -6  0   -9  15  0;...
    -3  6   3   -6  9   -21 6;...
    1   -3  2   0   -3  9   -6];
vL=(b*ones(1,6)).^[N+2,N+1,4,3,2,1];
vR=[N^3*(1-a),N^2*(1-a),N*(1-a),(1-a),N^2,N,1].';
prefac=1/(-6*(1-b)^4);
oa_b=vL*M1*vR*prefac;
ratio=abs(max(max(abs((vL.')*(vR.').*M1)))*prefac/oa_b);
if isnan(ratio)
    oa_b_ac=100;
else
    oa_b_ac=max(ratio*10^(-ndigits),N*(1-a));
end


 %[o_a_b_ac,o_ab_ac,oab_ac,oa_b_ac]
[accuracy,best]=min([o_a_b_ac,o_ab_ac,oab_ac,oa_b_ac,abz_ac]);
answs=[o_a_b,o_ab,oab,oa_b,abz];
bestAns=answs(best);

% if abs((brute-bestAns)/brute)>0.01
%     sprintf('a=%g, b=%g, Nn=%d, best=%d, bestAns=%g, brute=%g',a,b,N,best,bestAns,brute)
%     [o_a_b_ac,o_ab_ac,oab_ac,oa_b_ac,abz_ac]
%     error('innacuaracy found')
%     formatl long
%     a
%     b
% end

%out=[brute,o_a_b,o_ab,oab,oa_b,bestAns]; 
% if accuracy>0.1 && bestAns>1e-10
%     sprintf('a=%g, b=%g, best=%d, bestAns=%g',a,b,best,bestAns)
%     [o_a_b_ac,o_ab_ac,oab_ac,oa_b_ac,abz_ac]
%     error('innacuaracy found')    
% end
out=bestAns;

if isnan(out)
    error('NaN encountered')
end

end


