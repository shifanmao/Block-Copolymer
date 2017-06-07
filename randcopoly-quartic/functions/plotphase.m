clear;

s = load('solvent_spin');


% fprintf(outfile,'%.2f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', ...
%                  NM, LAM,  FA,   PHIP, CHIABS*NM,KSS*RM,EIGS(1)*NM,EIGS(2)*NM,CHIABKS*NM);

NM = 0.1;
FA = 0.3;
PHIP = 0.1;

IND = find(s(:,1) == NM & s(:,3) == FA & s(:,4)==PHIP);
LAMV = s(IND,2);
CHIABS = s(IND,5);
CHIABKS = s(IND,9);

figure;hold;set(gca,'fontsize',18)
plot(LAMV,CHIABS*PHIP,'k-','linewidth',2);
plot(LAMV(CHIABKS>0),CHIABKS(CHIABKS>0)*PHIP,'k--','linewidth',2);
xlabel('\lambda')
ylabel('4\chivf_Af_B\phi_P')
box on

% plot a straight line
CHIABV = [0:0.2:0.8];
for ii = 1:length(CHIABV)
    COL = (ii-1) / (length(CHIABV)-1);
    plot(-0.75,31.59*CHIABV(ii)*PHIP,'.','markersize',15,'color',[COL 0 1-COL])
end

% plot a straight line
CHIABV = [0:0.2:0.8];
for ii = 1:length(CHIABV)
    COL = (ii-1) / (length(CHIABV)-1);
    plot(0,20.0615*CHIABV(ii)*PHIP,'.','markersize',15,'color',[COL 0 1-COL])
end


%%%%%%%%%

s = load('solvent_spin');

NM = 0.1;
FA = 0.3;
PHIPV = linspace(0.1,1-1e-2,10);

figure;hold;set(gca,'fontsize',18)
for ii = 1:length(PHIPV)
    COL = (ii-1) / (length(PHIPV)-1);
    PHIP = PHIPV(ii);
    
    IND = find(s(:,1) == NM & s(:,3) == FA & abs(s(:,4)-PHIP)<1e-2);
    LAMV = s(IND,2);
    CHIABS = s(IND,5)
    CHIABKS = s(IND,9);

    plot(LAMV,CHIABS*PHIP,'k-','linewidth',2)
    %,'color',[0 0 COL]);
end









s = load('solvent_spin');
NM = 100;
IND = find(s(:,1) == NM & s(:,2) == 0);

PHIPV = s(IND,4);
FAV = s(IND,3);
CHIABS = s(IND,5);

PHIPV_unique = flipud(unique(PHIPV));
FAV_unique = unique(FAV);
Z = zeros(length(PHIPV_unique),length(FAV_unique));
for ii = 1:length(PHIPV_unique)
    for jj = 1:length(FAV_unique)
        index = find(PHIPV == PHIPV_unique(ii)  & FAV == FAV_unique(jj));
        Z(ii,jj) = CHIABS(index);
%         Z(ii,jj) = CHIABS(index)*NM*4*FAV_unique(jj)*(1-FAV_unique(jj));
    end
end

figure;imagesc(FAV_unique,PHIPV_unique,Z)
% figure;imagsc(FAV_unique,PHIPV_unique,Z,'EdgeColor','none')