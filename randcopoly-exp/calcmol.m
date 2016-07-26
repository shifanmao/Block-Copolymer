function [FA,rm]=calcmol(MW,WT)
% calculates mol percent (FA) and effective monomer length (rm)
% of ODPA-AP5F-PEG random copolymer system
% MW = Molecular Weight (g/mol) of PEG
% WT = Weight Percentage of PEG
if MW == 1500
    RB = 44.6;
elseif MW == 900
    RB = 39;
end

FA = 608*WT/(MW-(MW+306)*WT+608*WT);
rm = 30*(1-FA)+RB*FA;
end