% calculate lamellar spacing (D) in unit of 2lp
addpath('../functions')
addpath('../chainstats')
addpath('../misc')
addpath('../chainstats/eigcalc')
addpath('../chainstats/integrals')

%% Same N, different C simulations
clear;
% chain concentration
G = [20:8:52]; % total number of beads per chain

for N = logspace(1,4,4)
    NV = N*ones(1,length(G));

    % simulation parameters
    Dsim = 18;
    v = 0.1;

    % spacing in unit of 2lp
    [chis,ks,~]=spinodal(NV,0.5);
    D = 2*pi./ks;
    lptimes2 = Dsim./D;

    % end-to-end in unit of 2lp
    R2 = NV-(1/2)*(1-exp(-2*NV));
    R = sqrt(R2).*lptimes2;

    % concentration
    C = R.^3./(G*v);
    
    % renormalized spinodal
    [chit,phase]=spinodalRG(N,C,0.5);
    chit=reshape(chit,1,length(C));

    % discretization
    eps = NV./G;
    filename = sprintf('N%.2f',log10(N));
    dlmwrite(filename,[NV;log10(C);G/2;lptimes2;eps;chis.*NV;chit.*NV]','precision','%.3f')
end

%% Same C, different N simulations
clear;
% chain concentration
G = [20:8:52]; % total number of beads per chain

for C = logspace(1,4,4)
    CV = C*ones(1,length(G));

    % simulation parameters
    Dsim = 36;
    v = 0.1;

    % concentration
    R = power(CV.*G*v,1/3);
    fun = @(N) r2(N)-R.^2;
    N = lsqnonlin(fun,ones(1,length(R)),zeros(1,length(R)));

    % spacing in unit of 2lp
    [chis,ks,~]=spinodal(N,0.5);
    D = 2*pi./ks;
    lptimes2 = Dsim./D;
    
    % renormalized spinodal
    chit = zeros(1,length(N));
    for ii = 1:length(N)
        [chit(ii),phase]=spinodalRG(N(ii),C,0.5);
    end

    % discretization
    eps = N./G;
    filename = sprintf('C%.2f',log10(C));
    dlmwrite(filename,[N;log10(CV);G/2;lptimes2;eps;chis.*N;chit.*N]','precision','%.3f')
end