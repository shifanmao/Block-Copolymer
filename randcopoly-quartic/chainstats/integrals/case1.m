function S4=case1(N,NM,R1,R12,R3,F)   
% Case 1: J1==J2==J3==J4
% if separate return both valeq and valne
    E1=R1;
    E12=R12;
    E3=R3;
    
    % on same monomer
%     offset=10^-12;
%     if abs(E1)<offset
%         E1=offset;
%     end
%     if abs(E12)<offset
%         E12=offset;
%     end
%     if abs(E3)<offset
%         E3=offset;
%     end   
%{
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(NM*(E1+E12+E3)+5)/120;
    elseif max(abs([E1,E12]))<MIN
        valeq=chicken(E1,E12,E3,NM);
    elseif max(abs([E12,E3]))<MIN
        valeq=chicken(E12,E3,E1,NM);
    elseif max(abs([E1,E3]))<MIN
        valeq=chicken(E1,E3,E12,NM);
    elseif ((abs(E1-E12)<MIN) && (abs(E12-E3)<MIN))
        valeq=(1/2).*E1.^(-4).*((-2).*(3+E1.*NM)+exp(1).^(E1.*NM).*(6+E1.*NM.*(( ...
         -4)+E1.*NM)));
    elseif (abs(E1-E12)<MIN )%&& abs(E12-E3)>=MIN)
        valeq=E1.^(-3).*(E1+(-1).*E3).^(-2).*E3.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
          NM)).*E3.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E3.^2.*NM+E1.^3.*((-1)+ ...
          exp(1).^(E3.*NM)+(-1).*E3.*NM)+(-1).*E1.*E3.^2.*((-3)+E3.*NM+exp( ...
          1).^(E1.*NM).*(3+E3.*NM)));    
    elseif (abs(E1-E3)<MIN )%&& abs(E1-E12)>=MIN)
%         valeq=E1.^(-3).*(E1+(-1).*E12).^(-2).*E12.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
%           NM)).*E12.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E12.^2.*NM+E1.^3.*((-1)+ ...
%           exp(1).^(E12.*NM)+(-1).*E12.*NM)+(-1).*E1.*E12.^2.*((-3)+E12.*NM+ ...
%           exp(1).^(E1.*NM).*(3+E12.*NM)));
        valeq=((-3*E1*expl(4,NM*E1)+NM*E1^2*expl(3,NM*E1))*E12^2+...
          (2*expl(4,NM*E1)-NM*E1*expl(3,NM*E1))*E12^3+...
          E1^3*expl(4,NM*E12))/(E1^3*E12^2*(E1-E12)^2);
    elseif (abs(E12-E3)<MIN )%&& abs(E1-E3)>=MIN)
       valeq=E1.^(-2).*(E1+(-1).*E12).^(-2).*E12.^(-3).*(((-1)+exp(1).^(E1.*NM) ...
          ).*E12.^3+(-1).*E1.*E12.^3.*NM+E1.^2.*E12.*(3+2.*E12.*NM+exp(1).^( ...
          E12.*NM).*((-3)+E12.*NM))+(-1).*E1.^3.*(2+E12.*NM+exp(1).^(E12.* ...
          NM).*((-2)+E12.*NM)));
    else
       valeq=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*NM)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*NM+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*NM)).*E3.^2+E12.*E3.^2.*NM+E12.^2.*(( ...
          -1)+exp(1).^(E3.*NM)+(-1).*E3.*NM))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          NM)).*E3.^3+(-1).*E12.*E3.^3.*NM+E12.^3.*(1+(-1).*exp(1).^(E3.*NM) ...
          +E3.*NM)));
    end
%}
    valeq=case1Int(E1,E12,E3,NM);
    
    if valeq<0 && isreal(E1) && isreal(E12) && isreal(E3)
        sprintf('E1=%g, E12=%g, E3=%g',E1,E12,E3)
        error('shouldnt be negitive')
    end

    %sprintf('valeq=%g+%gi',real(valeq),imag(valeq))
    S4=zeros(2,2,2,2);
    S4(1,1,1,1)=S4(1,1,1,1)+F(1)*N*valeq;
    S4(2,2,2,2)=S4(2,2,2,2)+F(2)*N*valeq;

end
