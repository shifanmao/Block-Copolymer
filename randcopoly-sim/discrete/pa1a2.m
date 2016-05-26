function [PA1A2]=pa1a2(A1,A2,FA,LAM,J1,J2)

PAA=FA*(1-LAM)+LAM;
PBB=FA*(LAM-1)+1;
PBA=1-PAA;
PAB=1-PBB;
P=[PAA PAB;PBA PBB];
F=[FA 0;0 1-FA];

if A1==1
    FA1=FA;
else
    FA1=1-FA;
end

if A2==1
    FA2=FA;
else
    FA2=1-FA;
end

if J1 == J2
    if A1==A2
        PA1A2=FA1-FA1*FA2*0;
    else
        PA1A2=-FA1*FA2*0;
    end
else
    POW=abs(J1-J2);
    PMAT=(P^POW)*F;
    PA1A2=PMAT(A1,A2)-FA1*FA2*0;
end
