function heatOfReaction = shiftHeat(temperature)
%shiftHeat Calculate the standard heat of the water-gas shift reaction
%   temperature = temperature in K
%   heatOfReaction = standard heat of reaction in kJ/mol

% last modified on 2/14/18

    % species thermo
    [~, HsensCO, ~, dHformCO, ~] = thermoForCO(temperature);
    [~, HsensCO2, ~, dHformCO2, ~] = thermoForCO2(temperature);
    [~, HsensH2O, ~, dHformH2O, ~] = thermoForH2Ov(temperature);
    [~, HsensH2, ~, dHformH2, ~] = thermoForH2(temperature);
    
    % heat of reaction
    dHrxn298 = dHformCO2 + dHformH2 - dHformCO - dHformH2O;
    heatOfReaction = -HsensCO - HsensH2O + dHrxn298 + HsensCO2 + HsensH2;
end % of shiftHeat

