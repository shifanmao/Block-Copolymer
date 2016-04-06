function [val]=s3rr(N,G,LAM,FA,Q1,Q2,Q3)

val=zeros(2,2,2);
DISC=1;

if DISC==1

    for I1=1:N*G
        for I2=1:N*G
            for I3=1:N*G
              Qtot=norm(I1*Q1+I2*Q2+I3*Q3);
              pgc=zeros(2,2,2);
              
              for A1=1:2
                  for A2=1:2
                      for A3=1:2
                      if Qtot<1e-5
                          pgc(A1,A2,A3)=1;
                      else
                          pgc(A1,A2,A3)=sin(Qtot)/(Qtot);
                      end
                      val(A1,A2,A3)=val(A1,A2,A3)+pa1a2a3(A1,A2,A3,FA,LAM,...
                              floor((I1-1)/G),floor((I2-1)/G),floor((I3-1)/G))*pgc(A1,A2,A3);
                      end
                  end
              end
              
            end
        end
    end
else
    disp('no analytical Rigid Rod Solution')
end
end
 
function [PA1A2A3]=pa1a2a3(A1,A2,A3,FA,LAM,J1,J2,J3)
    F=[FA,1-FA];
    DF=[1,-1];
    [IND,ORD]=sort([J1,J2,J3]);
    A=[A1 A2 A3];
    OM=[A(ORD(1)),A(ORD(2)),A(ORD(3))];
    J1=IND(1);J2=IND(2);J3=IND(3);
    A1=OM(1);A2=OM(2);A3=OM(3);
    
    PA1A2A3=F(A1)*(F(A2)+DF(A1)*DF(A2)*(1-F(A1))*LAM^(J2-J1))*...
                  (F(A3)+DF(A2)*DF(A3)*(1-F(A2))*LAM^(J3-J2));
end

% function [PA1A2A3]=pa1a2a3(A1,A2,A3,FA,LAM,J1,J2,J3)
% 
% PAA=FA*(1-LAM)+LAM;
% PBB=FA*(LAM-1)+1;
% PBA=1-PAA;
% PAB=1-PBB;
% P=[PAA PAB;PBA PBB];
% F=[FA 0;0 1-FA];
% 
% [IND,ORD]=sort([J1,J2,J3]);
% A=[A1 A2 A3];
% OM=[A(ORD(1)),A(ORD(2)),A(ORD(3))];
% 
% if (IND(1) == IND(2)) && (IND(2)==IND(3))
%     if OM(1)==OM(2) && OM(2)==OM(3)
%         PA1A2A3=F(OM(1),OM(1));
%     else
%         PA1A2A3=0;
%     end
% elseif (IND(1) == IND(2)) && (IND(2)~=IND(3))
%     if OM(1)==OM(2)
%         POW2=abs(IND(3)-IND(2));
%         PMAT2=P^POW2*F;
%         PA1A2A3=PMAT2(OM(2),OM(3));
%     else
%         PA1A2A3=0;
%     end
% elseif (IND(1) ~= IND(2)) && (IND(2)==IND(3))
%     if OM(2)==OM(3)
%         POW1=abs(IND(2)-IND(1));
%         PMAT1=(P^POW1)*F;
%         PA1A2A3=PMAT1(OM(1),OM(2));
%     else
%         PA1A2A3=0;
%     end
% else
%     POW1=abs(IND(2)-IND(1));
%     POW2=abs(IND(3)-IND(2));
%     PMAT1=(P^POW1)*F;
%     PMAT2=P^POW2;
%     PA1A2A3=PMAT1(OM(1),OM(2))*PMAT2(OM(2),OM(3));
% end
% end
