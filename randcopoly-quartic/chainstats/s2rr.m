function s2=s2rr(N,G,LAM,FA,K)

s2=zeros(2,2);
DISC=1;

if DISC==1

    KJ=abs(K);
    for I1=1:N*G
        for I2=1:N*G
          sep = abs(I1-I2);
          pgc = zeros(2,2);
          if sep==0
              pgc(1,1)=1*pa1a2(1,1,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(1,2)=1*pa1a2(1,2,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(2,1)=1*pa1a2(2,1,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(2,2)=1*pa1a2(2,2,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
          else
              pgc(1,1)=sin(KJ*sep)/(KJ*sep)*pa1a2(1,1,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(1,2)=sin(KJ*sep)/(KJ*sep)*pa1a2(1,2,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(2,1)=sin(KJ*sep)/(KJ*sep)*pa1a2(2,1,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
              pgc(2,2)=sin(KJ*sep)/(KJ*sep)*pa1a2(2,2,FA,LAM,floor((I1-1)/G),floor((I2-1)/G));
          end
          s2(1,1)=s2(1,1)+pgc(1,1);
          s2(1,2)=s2(1,2)+pgc(1,2);
          s2(2,1)=s2(2,1)+pgc(2,1);
          s2(2,2)=s2(2,2)+pgc(2,2);
        end
    end
    
else

    disp('No analytical solution for continuous case');
    
end
end
 
function [PA1A2]=pa1a2(A1,A2,FA,LAM,J1,J2)
    F=[FA,1-FA];
    DF=[1,-1];
    [IND,ORD]=sort([J1,J2]);
    A=[A1 A2];
    OM=[A(ORD(1)),A(ORD(2))];
    J1=IND(1);J2=IND(2);
    A1=OM(1);A2=OM(2);
    
    PA1A2=F(A1)*(F(A2)+DF(A1)*DF(A2)*(1-F(A1))*LAM^(J2-J1));
end

% 
% function [PA1A2]=pa1a2(A1,A2,FA,LAM,J1,J2)
% 
% PAA=FA*(1-LAM)+LAM;
% PBB=FA*(LAM-1)+1;
% PBA=1-PAA;
% PAB=1-PBB;
% P=[PAA PAB;PBA PBB];
% F=[FA 0;0 1-FA];
% 
% if A1==1
%     FA1=FA;
% else
%     FA1=1-FA;
% end
% 
% if A2==1
%     FA2=FA;
% else
%     FA2=1-FA;
% end
% 
% if J1 == J2
%     if A1==A2
%         PA1A2=FA1;
%     else
%         PA1A2=0;
%     end
% else
%     POW=abs(J1-J2);
%     PMAT=(P^POW)*F;
%     PA1A2=PMAT(A1,A2)-FA1*FA2*0;
% end
