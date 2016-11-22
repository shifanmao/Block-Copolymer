function [REP,CHI,NTIME,NREP] = replica(chi,nodes)

% load data
[NTIME,NREP] = size(chi);
NREP = NREP-1;
TIME = 1:NTIME;

whos chi nodes
NODE = nodes(1:NTIME,2:end);
CHIV = chi(1:NTIME,2:end);
whos NODE CHIV

% order the replicas with nodes
REP = zeros(NTIME,NREP);  % 
CHI = zeros(NTIME,NREP);  % CHI
for ii = TIME
  for jj = 1:NREP
    TIND = find(TIME==ii);
    CIND = find(NODE(TIND,:) == jj);
    
    REP(ii,jj) = CIND;
    CHI(ii,jj) = CHIV(TIND,CIND);
  end
end
