clear;

filename = '../data/figureS4A';

N = 100;
NM = 1e-2;
FA = 0.6;

RM = sqrt(r2(NM));
PHIPV = [0.1, 0.9];

LAMMIN = max([-(1-FA)/FA, -FA/(1-FA)]);
LAMV = linspace(LAMMIN, 1, 51);

CHIABSV = zeros(length(PHIPV), length(LAMV));
KSV = zeros(length(PHIPV), length(LAMV));

for ii = 1:length(PHIPV)
    COL = (ii-1)/(length(PHIPV)-1);
    PHIP = PHIPV(ii);

    for jj = 1:length(LAMV)
        LAM = LAMV(jj)
        [CHIABSV(ii, jj), KSV(ii, jj), ~, ~] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
    end

%     figure(1);plot(LAMV, CHIABSV(ii, :)*4*FA*(1-FA)*PHIP*NM, ...
%         'color', [COL 0 1-COL]);
%     figure(2);plot(LAMV, RM*KSV(ii, :), ...
%         'color', [COL 0 1-COL])
end

% save(filename, 'CHIABSV', 'KSV', 'LAMV', 'PHIPV', 'NM');







% clear;
% 
% N = 100;
% NM = 0.01;
% 
% RM = sqrt(r2(NM));    
% FAV = linspace(0.1, 0.9, 51);
% % PHIPV = 0.1:0.1:0.9;
% PHIPV = [.1, .9];
% LAM = -0.75;
% FAV = linspace(-LAM/(1-LAM), 1/(1-LAM), 21);
% 
% CHIABSV = zeros(length(PHIPV), length(FAV));
% KSV = zeros(length(PHIPV), length(FAV));
% 
% for ii = 1:length(PHIPV)
%     COL = (ii-1)/(length(PHIPV)-1);
%     PHIP = PHIPV(ii)
% 
%     for jj = 1:length(FAV)
%         FA = FAV(jj)
%         [CHIABSV(ii, jj), KSV(ii, jj), ~, ~] = spinodal_wsolvent(N,NM,LAM,FA,PHIP);
%     end
% 
% %     figure(1);plot(LAMV, CHIABSV(ii, :)*4*FA*(1-FA)*PHIP*NM, ...
% %         'color', [COL 0 1-COL]);
% %     figure(2);plot(LAMV, RM*KSV(ii, :), ...
% %         'color', [COL 0 1-COL])
% end


% clear;
% filename1 = '../data/figureS5A';
% load(filename1);
% 
% figure;hold
% plot(FAV, CHIABSV(:, 11))
% plot(FAV, CHIABSV(:, 21))
% plot(FAV, CHIABSV(:, 31))
% plot(FAV, CHIABSV(:, 41))
% plot(FAV, CHIABSV(:, 51))