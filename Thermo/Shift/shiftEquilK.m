function equilK = shiftEquilK(temperature)
%shiftEquilK Calculate WGS equilibrium constant
%   temperature = temperature in K
%   equilK = equilibrium constant

% last modified on 4/5/18

    % gas constant
    idealR = 8.314E-3; % kJ/mol/K
    
    % species thermo
    [~, HsensCO, absSCO, dHformCO, ~] = thermoForCO(temperature);
    [~, HsensCO2, absSCO2, dHformCO2, ~] = thermoForCO2(temperature);
    [~, HsensH2O, absSH2O, dHformH2O, ~] = thermoForH2Ov(temperature);
    [~, HsensH2, absSH2, dHformH2, ~] = thermoForH2(temperature);
    
    % Enthalpy, entropy and free energy of reaction at 298K
    dHrxn298 = dHformCO2 + dHformH2 - dHformCO - dHformH2O;
    
    % Enthalpy and free energy of reaction at specified T
    dHrxnT = -HsensCO - HsensH2O + dHrxn298 + HsensCO2 + HsensH2;
    dSrxnT = absSCO2 + absSH2 - absSCO - absSH2O;
    dGrxnT = dHrxnT - temperature*dSrxnT/1000;
    
    % Equilibrium constant at specified T
    equilK = exp(-dGrxnT/idealR/temperature);
end

