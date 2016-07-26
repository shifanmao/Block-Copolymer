function S4=s4wlc(NM,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer)
%S4 is a four point correlation function
S4=0;
MIN=1e-5;

% Begin calculation of s4
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    
    disp(['sum(Q)=',num2str(sum(power(Q1+Q2+Q3+Q4,2)))])
    error('Wavevectors must add up to zero from translational invariance')

else

    % Evaluate the quantities for s4 calculation
    Q12=Q1+Q2;
    Q13=Q1+Q3;
    Q14=Q1+Q4;
    Q23=Q2+Q3;
    Q24=Q2+Q4;
    Q34=Q3+Q4;

    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));
    Q4MAG=sqrt(sum(power(Q4,2)));
    
    R1=Eigenvalues(norm(Q1MAG),ORDEig,ORDL);
    R2=Eigenvalues(norm(Q2MAG),ORDEig,ORDL);
    R3=Eigenvalues(norm(Q3MAG),ORDEig,ORDL);
    R4=Eigenvalues(norm(Q4MAG),ORDEig,ORDL);
    
    Q12MAG=sqrt(sum(power(Q12,2)));
    Q13MAG=sqrt(sum(power(Q13,2)));
    Q14MAG=sqrt(sum(power(Q14,2)));
    Q23MAG=sqrt(sum(power(Q23,2)));
    Q24MAG=sqrt(sum(power(Q24,2)));
    Q34MAG=sqrt(sum(power(Q34,2)));

    R12=Eigenvalues(Q12MAG,ORDEig,ORDL);
    R13=Eigenvalues(Q13MAG,ORDEig,ORDL);
    R14=Eigenvalues(Q14MAG,ORDEig,ORDL);
    R23=Eigenvalues(Q23MAG,ORDEig,ORDL);
    R24=Eigenvalues(Q24MAG,ORDEig,ORDL);
    R34=Eigenvalues(Q34MAG,ORDEig,ORDL);
    
    GL1=Residues(Q1MAG,R1,ORDEig,ORDL,NumLayer,1);
    GL2=Residues(Q2MAG,R2,ORDEig,ORDL,NumLayer,1);
    GL3=Residues(Q3MAG,R3,ORDEig,ORDL,NumLayer,1);
    GL4=Residues(Q4MAG,R4,ORDEig,ORDL,NumLayer,1);
    
    GLM12=Residues(Q12MAG,R12,ORDEig,ORDL,NumLayer,ORDL);
    GLM13=Residues(Q13MAG,R13,ORDEig,ORDL,NumLayer,ORDL);
    GLM14=Residues(Q14MAG,R14,ORDEig,ORDL,NumLayer,ORDL);
    GLM23=Residues(Q23MAG,R23,ORDEig,ORDL,NumLayer,ORDL);
    GLM24=Residues(Q24MAG,R24,ORDEig,ORDL,NumLayer,ORDL);
    GLM34=Residues(Q34MAG,R34,ORDEig,ORDL,NumLayer,ORDL);

    [YLM112,YLM443]=WignerD(Q1,Q12,Q4,ORDL);
    [~,YLM334]=WignerD(Q1,Q12,Q3,ORDL);
    [YLM113,YLM442]=WignerD(Q1,Q13,Q4,ORDL);
    [~,YLM224]=WignerD(Q1,Q13,Q2,ORDL);
    [YLM114,YLM332]=WignerD(Q1,Q14,Q3,ORDL);
    [~,YLM223]=WignerD(Q1,Q14,Q2,ORDL);
    [~,YLM441]=WignerD(Q2,Q23,Q4,ORDL);
    [~,YLM331]=WignerD(Q2,Q24,Q3,ORDL);
    [~,YLM221]=WignerD(Q3,Q34,Q2,ORDL);

    R1=NM*R1;
    R2=NM*R2;
    R3=NM*R3;
    R4=NM*R4;
    R12=NM*R12;
    R13=NM*R13;
    R14=NM*R14;
    R23=NM*R23;
    R24=NM*R24;
    R34=NM*R34;
    for M=0:(ORDL-1)
        if mod(ORDEig-M,2)==0
            Lmax=ORDEig-1;
        else
            Lmax=ORDEig-2;
        end

        for N1=0:Lmax
            for N2=M:Lmax
                for N3=0:Lmax
                    % PERMUTATIONS
                    % 1 12 4 112 443
                    % 1 12 3 112 334
                    % 1 13 4 113 442
                    % 1 13 2 113 224
                    % 1 14 3 114 332
                    % 1 14 2 114 223
                    % 2 23 4 223 441
                    % 2 23 1 223 114
                    % 2 24 3 224 331
                    % 2 24 1 224 113
                    % 3 34 2 334 221
                    % 3 34 1 334 112
    
                    % Case 1: A1=A2=A3=A4 (SAAAA,SBBBB)
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223);

                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113);

                    S4=S4_case1(S4,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221);
                    S4=S4_case1(S4,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112);
                end
            end
        end
    end
end
end

function S4=S4_case1(S4,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443)
% Case 1: A1=A2=A3=A4
[~,ORDL]=size(R12);
if M==0
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4=S4+2*S4_case1_int(1,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*...
                        YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                        *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1);
        end
    end
else
    for lam2=M:(ORDL-1)
        for lam3=M:(ORDL-1)
            S4=S4+4*S4_case1_int(1,R1(N1+1,1),R12(N2+1,M+1),R3(N3+1,1))*NM^4*real(...
                        YLM112(lam2+1,M+1)*conj(YLM443(lam3+1,M+1))...
                        *GL1(1,lam2+1,1,N1+1)*GLM12(lam2+1,lam3+1,M+1,N2+1)*GL3(lam3+1,1,1,N3+1));
        end
    end
end
end

% vvvvvvv helper functions vvvvvvvv
function valeq=S4_case1_int(NM,E1,E12,E3)
MIN=(10^-3)/NM;
vec=[E1,E12,E3];
nzeros=sum(abs(vec)<MIN);
if nzeros==3
    valeq=NM^4*(NM*(E1+E12+E3)+5)/120;
elseif nzeros==2
    [~,index]=max(abs(vec));
    others=vec; others(index)=[];
    valeq=chicken(others(1),others(2),vec(index),NM);
elseif nzeros==1
    [~,index]=min(abs(vec));
    others=vec; others(index)=[];
    if abs(others(1)-others(2))<MIN
        val=0.5*(others(1)+others(2));
        valeq=(-6*expl(3,NM*val)+2*NM*val*expl(2,NM*val))/(2*val^4);
    else
        valeq=(expl(4,NM*others(2))/(others(2)^3)...
            -expl(4,NM*others(1))/(others(1)^3))...
            /(others(2)-others(1));
    end
else
    dif=abs([E12-E3,E3-E1,E1-E12]);
    if max(dif)<MIN
        val=mean([E1,E12,E3]);
        valeq=(1/2)*val^(-4)*((-2)*(3+val*NM)+exp(val*NM)*(6+val*NM*(( ...
            -4)+val*NM)));
    elseif min(dif)>MIN
        valeq=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*NM)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*NM+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*NM)).*E3.^2+E12.*E3.^2.*NM+E12.^2.*(( ...
          -1)+exp(1).^(E3.*NM)+(-1).*E3.*NM))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          NM)).*E3.^3+(-1).*E12.*E3.^3.*NM+E12.^3.*(1+(-1).*exp(1).^(E3.*NM) ...
          +E3.*NM)));
    else
        [~,index]=min(dif);
        val=vec(index);
        others=vec; others(index)=[];
        oval=0.5*sum(others);
        valeq=((-3*oval*expl(4,NM*oval)+NM*oval^2*expl(3,NM*oval))*val^2+...
          (2*expl(4,NM*oval)-NM*oval*expl(3,NM*oval))*val^3+...
          oval^3*expl(4,NM*val))/(oval^3*val^2*(oval-val)^2);       
        
    end
    
end
end
function out=chicken(a,b,c,FA)
  out=(-FA^3*(FA*(a+b)+4)*c^4 ...
     -4*FA^2*(FA*(a+b)+3)*c^3 ...
    -12*FA^1*(FA*(a+b)+2)*c^2 ...
    +24*(expl(1,FA*c)-FA*(a+b))*c ...
    +24*expl(1,FA*c)*(a+b))/(24*c^5);
end