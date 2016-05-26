% clear
% 
% N=8;
% G=10;
% 
% FA=0.50;
% CHI=0.0/G;
% NK=100;
% 
% R2=G;   %if Gaussian chain
% KV = logspace(-2,2,NK)/sqrt(R2);
% val = zeros(NK,1);
% Sval = zeros(NK,1);
% for IK=1:NK
%     K=KV(IK);
%     LAM=-0.75;
% 
%     val(IK) = sk_gaussian(N,G,LAM,FA,K);
%     Sval(IK) = 1./(val(IK)-2*CHI);
% end
% 
% figure;
% semilogx(KV.*sqrt(R2),Sval)

clear

N=8;
G=10;
EPS=1.0;

FA=0.50;
CHI=0.0/G;
NK=100;

R2=G;   %if Gaussian chain
KV = logspace(-2,2,NK)/sqrt(R2);
val = zeros(NK,1);
Sval = zeros(NK,1);
for IK=1:NK
    K=KV(IK);
    LAM=-0.75;

    val(IK) = sk(N,G,LAM,FA,EPS,K);
    Sval(IK) = 1./(val(IK)-2*CHI);
end

figure;
loglog(KV.*sqrt(R2),Sval)