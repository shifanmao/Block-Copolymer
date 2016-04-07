function val=gamma4rep(N,NM,LAM,FA,k,Q1,Q2,Q3,Q4)
%% This function calculates the quartic order coefficient
% in free energy expansion of copolymer melts
% Usage: val=gamma4(N,FA,k,Q1,Q2,Q3,Q4)
% Return:
%    val, Quartic order expansion coeffiecient of free energy
% Parameters:
%    N, number of Kuhn steps
%    FA, fraction of A type monomers
%    k, magnitude of wavevector in unit of 1/contour length
%       k input can be a vector
%    Q1, Q2, Q3, Q4, wavevectors. Only relative orientations
%       matter, the magnitudes of Qs are all of Q
% Example:
%    N=1000;
%    k=sqrt(logspace(0,3,50)/N);
%    FA=0.5;
%    Q1(1:3)=[1,0,0];
%    Q2(1:3)=transpose(rotz(pi)*Q1(1:3)');
%    Q3(1:3)=-Q2;
%    Q4(1:3)=-Q1;
%    G=gamma4(k,Q1,Q2,Q3,Q4,CHAIN,NM,FA);
%    figure;hold;set(gca,'fontsize',15);leg=[];
%    plot(k.^2.*NM,G);
%    xlabel('(kR)^2');ylabel('\Gamma_4')
% Shifan Mao 06/10/15

MIN=1e-10;
if norm(Q1+Q2+Q3+Q4)>=MIN
    disp('ERROR :: Qs must add up to zero')
    return
else
    %result to return
    val=zeros(length(k),1);

    for j=1:length(k)
        %wave vectors
        Q1=Q1./norm(Q1)*k(j);
        Q2=Q2./norm(Q2)*k(j);
        Q3=Q3./norm(Q3)*k(j);
        Q4=Q4./norm(Q4)*k(j);

        % calculate single chain correlation functions
        if NM>=1e4  % Gaussian chain limit
            % Gaussian chain correlations
            s4 = s4gcrep(N,NM,LAM,FA,Q1,Q2,Q3,Q4);
        elseif NM<=1e-4  % Rigid rod limit
            % Rigid rod correlations
%             s4 = s4rrrep(N,NM,LAM,FA,Q1,Q2,Q3,Q4);
        else
            % Worm-like chain correlations
%             s4 = s4wlcrep(N,NM,LAM,FA,Q1,Q2,Q3,Q4);
        end
        s2inv=s2inverse(N,NM,LAM,FA,k(j));
        
        % use replica coupling term
        g4 = s4;

        M=combinator(2,4,'p','r');
        for I = 1:length(M)
            val(j) = val(j) + g4(M(I,1),M(I,2),M(I,3),M(I,4))*...
                    (s2inv(M(I,1),1)-s2inv(M(I,1),2))*...
                    (s2inv(M(I,2),1)-s2inv(M(I,2),2))*...
                    (s2inv(M(I,3),1)-s2inv(M(I,3),2))*...
                    (s2inv(M(I,4),1)-s2inv(M(I,4),2));
        end
    end
end

val=val*power(N*NM,3);
end