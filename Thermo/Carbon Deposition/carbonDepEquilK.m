function [KCH4, KCO, KH2O] = carbonDepEquilK(temperature)
%carbonDepEquilK Calculate equilibrium constants for C deposition
%   temperature = temperature in K
%   KCH4 = equilibrium constant for CH4 --> C + 2 H2
%   KCO = equilibrium constant for 2 CO --> C + CO2
%   KH2O = equilibrium constant for CO + H2 --> C + H2O

% last modified on 4/8/18

    % gas constant
    idealR = 8.314E-3; % kJ/mol/K
    
    % species thermo
    [~, HsensCO, absSCO, dHformCO, ~] = thermoForCO(temperature);
    [~, HsensCO2, absSCO2, dHformCO2, ~] = thermoForCO2(temperature);
    [~, HsensH2O, absSH2O, dHformH2O, ~] = thermoForH2Ov(temperature);
    [~, HsensH2, absSH2, dHformH2, ~] = thermoForH2(temperature);
    [~, HsensCH4, absSCH4, dHformCH4, ~] = thermoForCH4(temperature);
    
    % Enthalpies, entropies and free energies of reaction at 298K
    dHrxn298CH4 = 2*dHformH2 - dHformCH4;
    dHrxn298CO = dHformCO2 - 2*dHformCO;
    dHrxn298H2O = dHformH2O - dHformCO - dHformH2;
    
    % Enthalpies and free energies of reaction at specified T
    dHrxnTCH4 = 2*HsensH2 + dHrxn298CH4 - HsensCH4;
    dSrxnTCH4 = 2*absSH2 - absSCH4;
    dGrxnTCH4 = dHrxnTCH4 - temperature*dSrxnTCH4/1000;
    
    % Equilibrium constant at specified T
    equilKCH4 = exp(-dGrxnTCH4/idealR/temperature);
end

