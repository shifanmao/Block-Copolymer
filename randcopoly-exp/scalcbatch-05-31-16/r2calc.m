testnumv = [1,2,3,4,7,8,9,10,11,15,16,17,18];
cpnumv =   [1,1,1,2,2,1,1, 1, 1, 2, 1, 1, 2];

for I = 1:length(testnumv)
  r2calcfun(testnumv(I),cpnumv(I));
end
