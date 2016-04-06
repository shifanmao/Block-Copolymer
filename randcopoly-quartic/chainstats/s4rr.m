function [val]=s4rr(N,G,LAM,FA,Q1,Q2,Q3,Q4)

DISC=1;
val=zeros(2,2,2,4);

if DISC==1

    for I1=1:N*G
        for I2=1:N*G
            for I3=1:N*G
                for I4=1:N*G
                  Qtot=norm(I1*Q1+I2*Q2+I3*Q3+I4*Q4);
                  pgc=zeros(2,2,2,2);

              for A1=1:2
                  for A2=1:2
                      for A3=1:2
                          for A4=1:2
                          if Qtot<1e-5
                              pgc(A1,A2,A3,A4)=1;
                          else
                              pgc(A1,A2,A3,A4)=sin(Qtot)/(Qtot);
                          end
                          val(A1,A2,A3,A4)=val(A1,A2,A3,A4)+pgc(A1,A2,A3,A4)*pa1a2a3a4(A1,A2,A3,A4,FA,LAM,...
                                                            floor((I1-1)/G),...
                                                            floor((I2-1)/G),...
                                                            floor((I3-1)/G),...
                                                            floor((I4-1)/G));
                          end
                      end
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

function [PA1A2A3A4]=pa1a2a3a4(A1,A2,A3,A4,FA,LAM,J1,J2,J3,J4)
    F=[FA,1-FA];
    DF=[1,-1];
    [IND,ORD]=sort([J1,J2,J3,J4]);
    A=[A1 A2 A3 A4];
    OM=[A(ORD(1)),A(ORD(2)),A(ORD(3)),A(ORD(4))];
    
    J1=IND(1);J2=IND(2);J3=IND(3);J4=IND(4);
    A1=OM(1);A2=OM(2);A3=OM(3);A4=OM(4);
    
    PA1A2A3A4=F(A1)*(F(A2)+DF(A1)*DF(A2)*(1-F(A1))*LAM^(J2-J1))*...
                  (F(A3)+DF(A2)*DF(A3)*(1-F(A2))*LAM^(J3-J2))*...
                  (F(A4)+DF(A3)*DF(A4)*(1-F(A3))*LAM^(J4-J3));
end

% function [PA1A2A3A4]=pa1a2a3a4(A1,A2,A3,A4,FA,LAM,J1,J2,J3,J4)
% 
% PAA=FA*(1-LAM)+LAM;
% PBB=FA*(LAM-1)+1;
% PBA=1-PAA;
% PAB=1-PBB;
% P=[PAA PAB;PBA PBB];
% F=[FA 0;0 1-FA];
% 
% [IND,ORD]=sort([J1,J2,J3,J4]);
% A=[A1 A2 A3 A4];
% OM=[A(ORD(1)),A(ORD(2)),A(ORD(3)),A(ORD(4))];
% 
% POW1=abs(IND(2)-IND(1));
% POW2=abs(IND(3)-IND(2));
% POW3=abs(IND(4)-IND(3));
% PMAT1=(P^POW1)*F;
% PMAT2=P^POW2;
% PMAT3=P^POW3;
% PA1A2A3A4=PMAT1(OM(1),OM(2))*PMAT2(OM(2),OM(3))*PMAT3(OM(3),OM(4));
