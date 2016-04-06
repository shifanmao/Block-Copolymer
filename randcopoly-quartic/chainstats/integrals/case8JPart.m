function out=case8JPart(a,b,c,N)
%This function is used for quickly calculating:
%     out=0;
%     if 1
%         if N<31
%             for j4=4:N
%                 for j3=3:j4-1
%                     for j2=1:j3-1
%                         for j1=1:j2-1
%                             out=out+(c^(j4-j3))*(b^(j3-j2))*(a^(j2-j1));
%                         end
%                     end
%                 end
%             end
%         end
%         % return
%     end
%     brute=out;
% This functions accuarcy depends on it's inputs range.  I haven't worked
% hard on getting the relitive accuaracy for small output values.  Vary
% large N values can also cause problems (eg N~100,000).
% The variable "accuracy" can be used to estimate the accuaracy.
if abs(a*b*c)<1e-60
        out=0;
        return 
end
if N<4
    out=0;
    return
end

    goodEnough=0.001;
    
    if isreal(a) && isreal(b) && isreal(c)
        order=sort([c,b,a]);
        a=order(3);b=order(2);c=order(1); % a > b > c
    else
        %in the even that some of the values are complex
        % we use a slightly more complicated ordering process
        % see notes
        abc=[a,b,c];
        [junk,Index]=max(real(abc));
        a=abc(Index);
        abc(Index)=[];
        [junk,Index]=min(abs(abc-a));
        b=abc(Index);
        abc(Index)=[];
        c=abc;
        order=[c,b,a];
    end
    ndigits=14.5; %accuarcy of floating point numbers (all little conservative)
    
        % 1~a~b~c 
        oabc=-N*(N-1)*(N-2)*(N-3)*((N+1)*(3-a-b-c)-5)/120;
        oabc_ac=N*max(abs(1-[a,b,c]));
        if oabc_ac<goodEnough
            out=oabc;
            return
        end
        
        % 1~a~b, c~0
        oab_cz=c*(N-3)*(N-2)*(N-1)*(N*(a+b-2)+4)/24;
        %oab_cz=-c*(N*(b+a)+1-(N*(N+1)*((N^2-7*N+18)*(b+a-2)+4*(N-7)))/24)
        oab_cz_ac=max(abs([1-a,1-b,c]))*N;
        if oab_cz_ac<goodEnough
            out=oab_cz;
            return
        end       
        
        % 1~a b~c~0
        % oa_bcz=-c*b*(3*N-4*a+(N+1)*(18*a-3*N+N^2*a-7*N*a)/6-3);
        oa_bcz=b*c*((N+1)*(10*N+18*a+N^2*a-N^2-7*N*a-36)/6-4*a+10); 
        oa_bcz_ac=N*max(abs([1-a,b,c]));
        if oa_bcz_ac<goodEnough
            out=oa_bcz;
            return
        end  
        % a~b~c~z
        abcz=a*b*c*(N-3);
        abcz_ac=sum(abs([a,b,c]));
        if abcz_ac<goodEnough
            out=abcz;
            return
        end
        
        
        % 1~a !~ b~c
        M1=[...
            0   0   -6  18;...
            0   0   6   6;...
            0   0   0   0;...
            0   0   0   0;...
            -1  0   1   0;...
            3   -6  -9  0;...
            -3  12  -3  -18;...
            1   -6  11  -6;...
            0   0   0   0];
        M2=[...
            0   0   0   0;...
            0   -3  15  -18;...
            0   6   -12 -18;...
            0   -3   -3 0;...
            0   0   0   0;...
            0   0   0   0;...
            0   3   3   0;...
            0   -6  12  18;...
            0   3   -15 18];
        M3=[...
            0   0   6   -18;...
            0   0   -12 18;...
            0   0   6   0;...
            0   0   0   0;...
            0   3   -3  0;...
            0   -9  21  0;...
            0   9   -33 18;...
            0   -3  15  -18;...
            0   0   0   0];
        vL=(c*ones(1,9)).^[N+3,N+2,N+1,N,5,4,3,2,1];
        vR=[N^3;N^2;N;1];

        oa_bc=vL*(M1*(1-a)+ M2*(c-b) + M3)*vR/(6*(c-1)^5);
        
        if max(abs([N*(1-a), N*(b-c)/abs(b)]))< goodEnough && a-b>0.1
            out=oa_bc;
            return
        end        
        
        terms=(vL.')*(vR.').*(M1*(1-a)+ M2*(c-b) + M3)/(6*(c-1)^5);
        ratio=max(max(abs(terms)))/abs(oa_bc);
        if isnan(ratio)
            oa_bc_ac=100;
        else
            oa_bc_ac=max(abs([N*(1-a), N*(b-c)/abs(b), ratio*10^(-ndigits)]));
        end
        if oa_bc_ac<goodEnough
            out=oa_bc;
            return
        end       
        

        % 1 !~ a~b !~ c
         M1=[...
             0  0   N^2-5*N+6;...
             0  8*N-2*N^2-6   6*N-2*N^2;...
             N^2-3*N+2  4*N^2-8*N-8 N^2-N-2;...
             2*N-2*N^2+4    6-2*N^2 0;...
             N^2+N-2    0   0;...
             -2 -4  -4;...
             8  12  0;...
             -10    0   0];
         M2=[0  0   3   -6  3;...
             0  -2  3   0   -1;...
             0  0   0   0   0;...
             -1 0   0   0   1;...
             3  0   0   0   -3;...
             -3 0   -3  6   0;...
             1  2   -3  0   0];
         M3=[...
             0  0   -1  2   -1;...
             0  1   -1  -1  1;...
             0  -1  2   -1  0;...
             0  0   0   1   -1;...
             0  0   -2  1   1;...
             0  1   1   -2  0;...
             0  -1  1   0   0];

         vL=((b*ones(1,8)).^[N+4,N+3,N+2,N+1,N,5,4,3]);
         vR=[c^2;c;1];
         T1=c*(b-1)^(-2)*vL*M1*vR;
         terms1=max(max(abs((vL.')*(vR.').*M1)*c))/((1-b)^2);
         
         T2=-2*b^2*(c^3-c^N)*(b-2*c+b*c)/(c-1);
         terms12=max(terms1,abs(T2));

         T3=-c^2*(1-c)*(b^4-b^(N+1))*(1-b)^(-1)...
            +b^2*(1-b)*(c^4-c^(N+1))*(1-c)^(-1)...
            +b^2*c^2*(b-c)*(N-3);
         terms3=max(abs([-c^2*(1-c)*(b^4-b^(N+1))*(1-b)^(-1), ...
                           b^2*(1-b)*(c^4-c^(N+1))*(1-c)^(-1), ...
                           b^2*c^2*(b-c)*(N-3)]));
        
         vL=(b*ones(1,7)).^[N+2,N+1,N,4,3,2,1];
         vR=(c*ones(5,1)).^[N;3;2;1;0];
         T4=b*c*( (b-c)^2*(b-1)^3*(c-1)^2 )^(-1)*...
             vL*(M2+M3*N)*vR;
         terms4=b*c*( (b-c)^2*(b-1)^3*(c-1)^2 )^(-1)*...
                max(max(abs((vL.')*(vR.').*(M2+M3*N))));

        o_ab_c=(a-b)*( (T1-T2)/(2*(b-1)^2*(b-c)^3) + T3/(b*c*(b-c)*(b-1)^3*(c-1)))-T4;
        
        bigmax=max(abs([(a-b)*terms12/(2*(b-1)^2*(b-c)^3), ...
                        (a-b)*terms3/(b*c*(b-c)*(b-1)^3*(c-1)), ...
                        terms4]));
        ratio=bigmax/abs(o_ab_c);
        if isnan(ratio)
            o_ab_c_ac=100;
        else
            o_ab_c_ac=max(abs([N*(a-b)/abs(b), ratio*10^(-ndigits)]));
        end
        if o_ab_c_ac<goodEnough
            out=o_ab_c;
            return
        end 
        
        % 1 !~ a !~ b~c
        %This is the same as above with a and c switched
        a=order(1);b=order(2);c=order(3);
 
         M1=[...
             0  0   N^2-5*N+6;...
             0  8*N-2*N^2-6   6*N-2*N^2;...
             N^2-3*N+2  4*N^2-8*N-8 N^2-N-2;...
             2*N-2*N^2+4    6-2*N^2 0;...
             N^2+N-2    0   0;...
             -2 -4  -4;...
             8  12  0;...
             -10    0   0];
         M2=[0  0   3   -6  3;...
             0  -2  3   0   -1;...
             0  0   0   0   0;...
             -1 0   0   0   1;...
             3  0   0   0   -3;...
             -3 0   -3  6   0;...
             1  2   -3  0   0];
         M3=[...
             0  0   -1  2   -1;...
             0  1   -1  -1  1;...
             0  -1  2   -1  0;...
             0  0   0   1   -1;...
             0  0   -2  1   1;...
             0  1   1   -2  0;...
             0  -1  1   0   0];

         vL=((b*ones(1,8)).^[N+4,N+3,N+2,N+1,N,5,4,3]);
         vR=[c^2;c;1];
         T1=c*(b-1)^(-2)*vL*M1*vR;
         terms1=max(max(abs((vL.')*(vR.').*M1)*c))/((1-b)^2);
         
         T2=-2*b^2*(c^3-c^N)*(b-2*c+b*c)/(c-1);
         terms12=max(terms1,abs(T2));

         T3=-c^2*(1-c)*(b^4-b^(N+1))*(1-b)^(-1)...
            +b^2*(1-b)*(c^4-c^(N+1))*(1-c)^(-1)...
            +b^2*c^2*(b-c)*(N-3);
         terms3=max(abs([-c^2*(1-c)*(b^4-b^(N+1))*(1-b)^(-1), ...
                           b^2*(1-b)*(c^4-c^(N+1))*(1-c)^(-1), ...
                           b^2*c^2*(b-c)*(N-3)]));
        
         vL=(b*ones(1,7)).^[N+2,N+1,N,4,3,2,1];
         vR=(c*ones(5,1)).^[N;3;2;1;0];
         T4=b*c*( (b-c)^2*(b-1)^3*(c-1)^2 )^(-1)*...
             vL*(M2+M3*N)*vR;
         terms4=b*c*( (b-c)^2*(b-1)^3*(c-1)^2 )^(-1)*...
                max(max(abs((vL.')*(vR.').*(M2+M3*N))));

        o_a_bc=(a-b)*( (T1-T2)/(2*(b-1)^2*(b-c)^3) + T3/(b*c*(b-c)*(b-1)^3*(c-1)))-T4;
        
        bigmax=max(abs([(a-b)*terms12/(2*(b-1)^2*(b-c)^3), ...
                        (a-b)*terms3/(b*c*(b-c)*(b-1)^3*(c-1)), ...
                        terms4]));
        ratio=bigmax/abs(o_a_bc);
        if isnan(ratio)
            o_a_bc_ac=100;
        else
            o_a_bc_ac=max(abs([N*(a-b)/abs(b), ratio*10^(-ndigits)])); 
        end
        a=order(3);b=order(2);c=order(1); % put back in origional order
        if o_a_bc_ac<goodEnough
            out=o_a_bc;
            return
        end         
        
        % a~b~c !~ 1
        % I have measured w.r.t. c as this gives best results
        M1=[...
            0   0   0   0;...
            -1  6   -11 6;...
            3   -12 3   18;...
            -3  6   9   0;...
            1   0   -1  0;...
            0   0   0   0;...
            0   0   0   0;...
            0   0   -6  -6;...
            0   0   6   -18];

        M2=[...
            0   3   -15 18;...
            0   -9  33  -18;...
            0   9   -21 0;...
            0   -3  3   0;...
            0   0   0   0;...
            0   0   -6  0;...
            0   0   12  -18;...
            0   0   -6  18;...
            0   0   0   0];

        vL=(c*ones(1,9)).^[N+4,N+3,N+2,N+1,N,5,4,3,2];
        vR=[N^3;N^2;N;1];

        o_abc=vL*(M1*((b-a)+2*(c-b))+M2)*vR/(6*(c-1)^5);
        
        terms=(vL.')*(vR.').*(M1*((b-a)+2*(c-b))+M2)/(6*(c-1)^5);
        ratio = max(max(abs(terms)))/abs(o_abc);
        if isnan(ratio)
            o_abc_ac=100;
        else
            o_abc_ac=max(abs([N*(a-b)/abs(a), N*(b-c)/abs(a), ratio*10^(-ndigits)]));
        end
        if o_abc_ac<goodEnough
            out=o_abc;
            return
        end
        
        %  1~a~b !~ c
        M1=...
        [1,-2,-1,2,0;...
         -4,12,4,-12,0;...
         6,-24,6,36,0;...
         -4,20,-20,-20,24;...
         1,-6,11,-6,0];

        M2=...
        [0,2,-6,4,0;...
         0,-8,30,-22,0;...
         0,12,-54,54,0;...
         0,-8,42,-58,12;...
         0,2,-12,22,-12];

        vL=[c^5,c^4,c^3,c^2,c];
        vR=[N^4;N^3;N^2;N;1];

        oab_c=(vL*( ((a-b)+2*(1-a))*M1-2*M2 )*vR...
            -((a-b)+2*(1-a))*24*c^(N+2)...
            -24*c^(N+1)+24*c^(N+2) )...
            /(24*(c-1)^5);
        
        terms=(vL.')*(vR.').*( ((a-b)+2*(1-a))*M1-2*M2 )/(24*(c-1)^5);
        ratio = max(max(abs(terms)))/abs(o_abc);
        if isnan(ratio)
            oab_c_ac=100;
        else
            oab_c_ac = max(abs([N*(1-b), ratio*10^(-ndigits)]));
        end
        if oab_c_ac<goodEnough
            out=oab_c;
            return
        end        
        
        % 1~a !~ b !~ c
        M1=[...
            1   0   -1  -6;...
            -3  3   6   18;...
            3   -6  -3  -12;...
            -1  3   -2  6];
        M2=[...
            0   -3  3   6;...
            0   9   -15 -18;...
            0   -9  21  12;...
            0   3   -9  0];
        M3=[...
            -2  3   5   -6;...
            6   -15 -15 24;...
            -6  21  3   -36;...
            2   -9  7   6];
        M4=[...
            0   6   -12 0;...
            0   -18 48  -6;...
            0   18  -60 24;...
            0   -6  24  -18];
        M5=[...
            1   -3  2   0;...
            -3  12  -9  0;...
            3   -15 18  0;...
            -1  6   -11 6];
        M6=[...
            0   -3  9   -6;...
            0   9   -33 24;...
            0   -9  39  -36;...
            0   3   -15 18];
        M7=[...
            -18 18  -6;...
            18  -6  0;...
            -6  0   0];
        M8=[...
            12  -18 6;...
            -18 24  -6;...
            6   -6  0];
        vL=((c*ones(1,4)).^[5,4,3,2]);
        vR=((N*ones(4,1)).^[3;2;1;0]);
        MM=M1*(1-a)*b^3+...
           M2*b^3+...
           M3*(1-a)*b^2+...
           M4*b^2+...
           M5*(1-a)*b+...
           M6*b;
        alpha1=vL*MM*vR;
        alpha2=[b^3,b^2,b]*(M7*(1-a)+M8)*[c^(N+2);c^(N+1);c^N];

        alpha3=b*(b*(1-a)-b+1) * ((c^2-c^3)*b^N+(c^3-c^N)*b^2+(c^N-c^2)*b^3)/...
               ( c*(b-c)*(b-1)*(c-1) );
        oa_b_c=-(b-1)^(-3)*((alpha1+alpha2)/(6*c*(c-1)^4) -alpha3);  

        terms=(vL.')*(vR.')/(6*c*(b-1)^3*(c-1)^4);
        terms2=max(abs([(c^2-c^3)*b^N, (c^N-c^2)*b^3, (c^3-c^N)*b^2]))...
               *b/(c*(b-c)*(b-1)^4*(c-1));
        ratio=max(max(max(abs(terms))),abs(terms2))/abs(oa_b_c);
        if isnan(ratio)
            oa_b_c_ac=100;
        else
            oa_b_c_ac=max(abs([N*(1-a), ratio*10^(-ndigits)]));
        end
        if oa_b_c_ac<goodEnough
            out=oa_b_c;
            return
        end 
        
        % a !~ b~c~z  also includes 1~a !~ b~c~z
        a_bcz_ac=abs(b)+abs(c);
        a_bcz=a*b*c*((1-a)^(-1))*((-a^3+a^N)/(a^2*(1-a))+N-3);
        if abs((N-3)/((1-a)*a_bcz*10^(ndigits)))>a_bcz_ac && a > 0.999
            a_bcz=b*c*0.5*(4-2*N+(N-1)*(N-2)*a^(N-3));
        end
        
        
        % a~b~z !~ c
        abz_c_ac=abs(a)+abs(b);
        abz_c=a*b*c*((1-c)^(-1))*((-c^3+c^N)/(c^2*(1-c))+N-3);
        
        oa_bz_c=-b*(6*c+8*N*c^2+N^2*c-3*N*c^3-6*c^2+2*c^3-2*N^2*c^2+N^2*c^3-5*N*c-2*c^N)/(2*(c-1)^3);
        oa_bz_c_ac=max(abs([N^3*(1-a),b,10^(-ndigits)/((c-1)^3)]));
        
        
        % compair accuaracy of various approximations
        cut=0.001/(abs(a-b)*abs(b-c));
        accuracies=[oabc_ac,...
                              o_abc_ac,...
                              oab_c_ac,...
                              oa_bc_ac,...
                              oa_b_c_ac,...
                              o_ab_c_ac,...
                              o_a_bc_ac,...
                              abcz_ac,...
                              a_bcz_ac,...
                              abz_c_ac,...
                              cut,...   Calculate this one later for speed
                              oab_cz_ac,...
                              oa_bcz_ac,...
                              oa_bz_c_ac];
        [accuracy, best]=min(accuracies);
        if best==11
            % a !~ b !~ c
            if( min(abs([1-a,1-b,1-c]))>0.02)  %This option is for speede
                T1=(-b^2*(1-a)/((b-a)*(c-b))+a^2*(1-b)/((b-a)*(c-a))+1/(c-1))*(c^3-c^N)/(-c^2);
                T2=(a/(b-a)+1/(1-b))*(c-1)*(b^3-b^N)/(b*(c-b));
                T3=(1-b)*(1-c)*(a^3-a^N)/((b-a)*(c-a)*(1-a))-(N-3);
                a_b_c=-a*b*c*(T1+T2+T3)/((1-a)*(1-b)*(1-c));           
            else
            x1=log(a);
            x2=log(b);
            x3=log(c);

            outs=[0;0;0;0];
            v1=[x1,x1,x3,x3,x1-x3,3*x2-N*x2];
            v2=[expl(2,x1),expl(2,x1),expl(2,x3),expl(2,x3),...
                expl(2,x1)-expl(2,x3),expl(2,3*x2)-expl(2,N*x2)];
            outs(1)= hooperHumperdink(v1,v2);

            v1=[x2,x2,x3,x3,x2-x3,3*x1-N*x1];
            v2=[expl(2,x2),expl(2,x2),expl(2,x3),expl(2,x3),...
                expl(2,x2)-expl(2,x3),expl(2,3*x1)-expl(2,N*x1)];
            outs(2)=(-1)*hooperHumperdink(v1,v2);

            v1=[x1,x1,x2,x2,x1-x2,3*x3-N*x3];
            v2=[expl(2,x1),expl(2,x1),expl(2,x2),expl(2,x2),...
                expl(2,x1)-expl(2,x2),expl(2,3*x3)-expl(2,N*x3)];
            outs(3)=(-1)*hooperHumperdink(v1,v2);

            v1=[x1,x2,x3,x1-x2,x1-x3,x2-x3];
            v2=[expl(2,x1),expl(2,x2),expl(2,x3),...
                expl(2,x1)-expl(2,x2),...
                expl(2,x1)-expl(2,x3),...
                expl(2,x2)-expl(2,x3)];
            outs(4)=hooperHumperdink(v1,v2)*(3-N);
            out=sum(outs);

            a_b_c=out*a*b*c*( (a-1)^2*(b-1)^2*(c-1)^2*(a-b)*(a-c)*(b-c) )^(-1);
            end
        else
            a_b_c=0;
        end

        answs=[oabc,...
               o_abc,...
               oab_c,...
               oa_bc,...
               oa_b_c,...
               o_ab_c,...
               o_a_bc,...
               abcz,...
               a_bcz,...
               abz_c,...
               a_b_c,...
               oab_cz...
               oa_bcz,...
               oa_bz_c];
%         if 0.01<abs((answs(best)-brute)/brute)
%             sprintf('a=%g, b=%g, c=%g, 1-a=%g, a-b=%g, b-c=%g',a,b,c,1-a,a-b,b-c)
%             sprintf('brute=%g, best=%g',brute,answs(best))
%             sprintf('best #=%d',best)
%             sprintf('accuracy=%g',accuracy)
%             error('found innacuarcy')
%         end
            
        out=answs(best);
        if isnan(out)
            sprintf('a=%g+%gi,b=%g+%gi,c=%g+%gi,N=%d',...
                     real(a),imag(a),real(b),imag(b),real(c),imag(c),N)
            sprintf('element %d of',best)
            disp(answs)
            disp('accuracies:')
            disp(accuracies)
            error('incountered NaN in case8Jpart')
        end
end
function out=nearOneAndZero(a,b,c,N)
abc=[a,b,c];
rounded=real(round(abc));
accuracy=sum(abs(rounded-abc))*N^sum(rounded);

end
