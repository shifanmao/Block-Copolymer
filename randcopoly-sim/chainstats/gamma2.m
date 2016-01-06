function val=gamma2(M,NM,LAM,FA,Q,CHI)
% Function : gamma2
% This function calculates the quadratic order coefficient
% in free energy expansion of copolymer melts
% Usage: val=gamma2(M,NM,LAM,FA,Q,CHI)
% Return:
%     val, Quadratic order expansion coeffiecient of free energy
% of random copolymer melt
% Parameters:
%     M, number of monomers
%     NM, number of Kuhn steps per monomer
%     LAM, chemical correlation of random copolymers
%     FA, fraction of A type monomers
%     Q, magnitude of wavevector in unit of 1/contour length
%        Q input can be a vector
%     CHI, chemical incompatibility between A and B monomers, non-
%         dimensionalized by monomer volume v (CHI*v)
% Example:
%     k=logspace(-6,2,50);
%     g=gamma2(100,1e5,-.75,0.5,k,0.00005);
%     figure;loglog(k,1./g);
%     xlabel('k');ylabel('1/\Gamma_2')
% This function calculates the quadratic coefficient
% in free energy expansion of Gaussian chain, Wormlike Chain
% and Rigid rod
%
% Shifan Mao 06/10/15

if NM>=1e4  % Gaussian chain limit
    val = -2*CHI+NM*gamq2gc(M,NM,LAM,Q)/(FA*(1-FA));
elseif NM<=1e-3  % Rigid rod limit
else  % Worm-like chain
    val = -2*CHI+NM*gamq2wlc(M,NM,LAM,Q)/(FA*(1-FA));
end
end