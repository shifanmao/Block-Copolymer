function [gam3,gam4]=calcgamma(N,FAV,NQ)
%% calcgamma.m :: This code calculates the coefficients in free energy expansion of block
% copolymers. Coefficients are evaluated by finding the cubic and quartic
% order vertex functions at dominant peak given in quadratic order fluctuations
% Usage: [gam3,gam4]=calcgamma(N,FAV,NQ)
% Parameters:
%   CHAIN, type of polymer (CHAIN=1, Guassian; =2, WLC; =3 Rigid Rod)
%   NM, number of Kuhn steps of total chain
%   FAV, range of A-type monomer fractions
%   ORDEig, number of eigenvalues
%   ORDL, number of spherical harmonics
%   NumLayer, Number of residual layers
% Return:
%   chis, Flory-Huggins parameter at spinodal
%   ks, critical wavelength of quadratic density fluctuation
%   gam3, cubic order vertex constant
%   gam4, quartic order vertex function

filename=sprintf('data/gam%.2e.mat',N);
if exist(filename,'file')
    fprintf('Step 2: Loading vertices at N=%.2e\n',N)
    load(filename);
    gam3=GAM3;
    gam4=GAM4;
else
    % wavevectors for Gamma4 calculations
    Q1=zeros(3,NQ);
    Q2=zeros(3,NQ);
    Q3=zeros(3,NQ);
    Q4=zeros(3,NQ);

    if NQ==1
        % case1 :: theta=pi
        Q1(1:3,1)=[1,0,0];
        Q2(1:3,1)=rotz(pi)*Q1(1:3,1);
        Q3(1:3,1)=-Q2(1:3,1);
        Q4(1:3,1)=-Q1(1:3,1);

    elseif NQ==4
        % case1 :: theta=pi
        Q1(1:3,1)=[1,0,0];
        Q2(1:3,1)=rotz(pi)*Q1(1:3,1);
        Q3(1:3,1)=-Q2(1:3,1);
        Q4(1:3,1)=-Q1(1:3,1);

        % case2 :: theta=pi/3
        Q1(1:3,2)=[1,0,0];
        Q2(1:3,2)=rotz(pi/3)*Q1(1:3,2);
        Q3(1:3,2)=-Q1(1:3,2);
        Q4(1:3,2)=-Q2(1:3,2);

        % case3 :: theta=pi/2
        Q1(1:3,3)=[1,0,0];
        Q2(1:3,3)=rotz(pi/2)*Q1(1:3,3);
        Q3(1:3,3)=-Q1(1:3,3);
        Q4(1:3,3)=-Q2(1:3,3);

        % case4 :: gam4(1,2)
        Q1(1:3,4)=[-1,0,1];
        Q2(1:3,4)=[-1,0,-1];
        Q3(1:3,4)=[1,1,0];
        Q4(1:3,4)=[1,-1,0];
    else
        Qrot=linspace(0,pi/2,NQ);
        for IQ=1:NQ
            Qvec=Qrot(IQ);
            Q1(1:3,IQ)=[1,0,0];
            Q2(1:3,IQ)=rotz(Qvec)*Q1(1:3,IQ);
            Q3(1:3,IQ)=-Q1(1:3,IQ);
            Q4(1:3,IQ)=-Q2(1:3,IQ);
        end
    end

    % results to return :: cubic and quartic order coefficients
    gam3=zeros(length(FAV),1);
    gam4=zeros(length(FAV),NQ);

    % calculate spinodal and critical wavelength
    [~,ks]=spinodal(N,FAV);

    for IFA=1:length(FAV)
        FA=FAV(IFA);
        fprintf('Step 2: Calculating vertices at FA=%.2f, N=%.2e\n',FA,N)

        % calculate free energy coefficients
        gam3(IFA) = gamma3(N,FA,ks(IFA));

        for IQ=1:NQ
            K1=Q1(:,IQ);
            K2=Q2(:,IQ);
            K3=Q3(:,IQ);
            K4=Q4(:,IQ);
            gam4(IFA,IQ) = gamma4(N,FA,ks(IFA),K1,K2,K3,K4);
        end
    end
end
end