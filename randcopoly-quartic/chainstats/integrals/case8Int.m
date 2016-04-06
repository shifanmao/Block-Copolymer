function valeq=case8Int(E1,E12,E3,NM)

    % on same monomer
    MIN=(10^-4)/NM;    
    
    vec=[E1,E12,E3];
    nzeros=sum(abs(vec)<MIN);
    if nzeros==3
        valeq=NM^4;
    elseif nzeros==2
        if max(abs([E1,E12]))<MIN
            valeq=(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))...
                -2*NM*E12*E3*sinh(NM*E3));
        elseif max(abs([E12,E3]))<MIN
            valeq=(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))...
                -2*NM*E12*E1*sinh(NM*E1));
        elseif max(abs([E1,E3]))<MIN
            valeq=(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2 ...
                -NM*(E1+E3)*E12*sinh(NM*E12));
        else
            error('not an option')
        end
    elseif nzeros==1
        if abs(E1-E12)<MIN && abs(E3)<MIN
            valeq=4*NM^2*sinh(0.5*NM*E1)^2/(E1^2);
        elseif abs(E12-E3)<MIN && abs(E1)<MIN
            valeq=4*NM^2*sinh(0.5*NM*E3)^2/(E3^2);
        elseif abs(E1-E3)<MIN && abs(E12)<MIN
            valeq=16*sinh(0.5*NM*E1)^4/(E1^4);
        elseif abs(E1)<MIN
            valeq=-NM*(expl(1,NM*E12)-expl(1,NM*E3))*expl(1,-NM*E12)*expl(1,-NM*E3)/(-E12*E3*(E12-E3));
        elseif abs(E12)<MIN   
            valeq=16*sinh(0.5*NM*E1)^2*sinh(0.5*NM*E3)^2/(E1^2*E3^2);
        elseif abs(E3)<MIN
            valeq=-NM*(expl(1,NM*E12)-expl(1,NM*E1))*expl(1,-NM*E12)*expl(1,-NM*E1)/(-E12*E1*(E12-E1));
        else
            error('not an option')
        end
    else
        if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
           valeq=2.*E1.^(-2).*NM.^2.*((-1)+cosh(E1.*NM));
        elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
           valeq=exp(1).^((-1).*(E1+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
              E1.*NM)+(-1).*exp(1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1) ...
              .*(E1+(-1).*E3).^(-1).*E3.^(-1).*NM;
        elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
           valeq=16.*E1.^(-2).*(E1+(-1).*E12).^(-2).*sinh((1/2).*E1.*NM).^2.*sinh(( ...
            1/2).*(E1+(-1).*E12).*NM).^2;
        elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
           valeq=exp(1).^((-1).*(E1+E12).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
              E1.*NM)+(-1).*exp(1).^(E12.*NM)).*((-1)+exp(1).^(E12.*NM)).*E1.^( ...
              -1).*(E1+(-1).*E12).^(-1).*E12.^(-1).*NM;
        else
    %        valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
    %           .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
    %           1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
    %           .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);    
            valeq=2*( -coshl(4,NM*E1)+coshl(4,NM*E12)-coshl(4,NM*E3)-coshl(4,NM*(E1-E12))+...
               coshl(4,NM*(E1-E3))-coshl(4,NM*(E12-E3))+coshl(4,NM*(-E1+E12-E3)) )...
               *( (E12-E3)*E3*E1*(E12-E1) )^(-1);    
        end
    end 



%{
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
       valeq=NM^4;
    elseif max(abs([E1,E12]))<MIN
        valeq=(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))...
            -2*NM*E12*E3*sinh(NM*E3));
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))...
            -2*NM*E12*E1*sinh(NM*E1));
    elseif max(abs([E1,E3]))<MIN
        valeq=(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2 ...
            -NM*(E1+E3)*E12*sinh(NM*E12));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=2.*E1.^(-2).*NM.^2.*((-1)+cosh(E1.*NM));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1) ...
          .*(E1+(-1).*E3).^(-1).*E3.^(-1).*NM;
    elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
       valeq=16.*E1.^(-2).*(E1+(-1).*E12).^(-2).*sinh((1/2).*E1.*NM).^2.*sinh(( ...
        1/2).*(E1+(-1).*E12).*NM).^2;
    elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E12).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E12.*NM)).*((-1)+exp(1).^(E12.*NM)).*E1.^( ...
          -1).*(E1+(-1).*E12).^(-1).*E12.^(-1).*NM;
    else
%        valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
%           .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
%           1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
%           .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);    
        valeq=2*( -coshl(4,NM*E1)+coshl(4,NM*E12)-coshl(4,NM*E3)-coshl(4,NM*(E1-E12))+...
           coshl(4,NM*(E1-E3))-coshl(4,NM*(E12-E3))+coshl(4,NM*(-E1+E12-E3)) )...
           *( (E12-E3)*E3*E1*(E12-E1) )^(-1);    
    end
%}
end