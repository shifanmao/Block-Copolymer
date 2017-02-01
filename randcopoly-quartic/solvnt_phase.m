clear;

N = 100;
NM = 0.1;
FA = 0.5;
FP = 0.5;   % Fraction of polymers
LAM = -0.75;

% WAVEVECTORS
NK = 100;
RM = (r2(NM))^0.5; % Normalization factor
KV = logspace(-2,2,NK)/RM;  % Wavevector

% CHI PARAMETERS
CHIAB = 0/NM; % Flory-Huggins factor between A and B
CHIBS = 0/NM; % Flory-Huggins factor between B and S
CHIAS = 0;
CHIBA = CHIAB;
CHI = [CHIAS,CHIAB;CHIBA,CHIBS];
[EIG1,EIG2,EIGV1,EIGV2,THETA1,THETA2,KS1,KS2]=gamma2_solvent(N,NM,LAM,FA,KV,CHI,FP);

t = linspace(0,2*pi,100);
x = sin(t);
y = cos(t);

f = figure('position', [0, 0, 800, 800]);hold
plot(x,y,'k-')
plot([0, -EIGV1(1,1)], [0, -EIGV1(1,2)],'r-','linewidth',2)
plot([0, -EIGV1(1,1)], [0, -EIGV1(1,2)],'-','linewidth',2,'color',rand(1,3))
box on