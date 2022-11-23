function equilK = steamReformingEquilK(temperature)
%shiftEquilK Calculate steam reforming equilibrium constant
%   CH4 + H2O <--> CO + 3 H2
%   temperature = temperature in K
%   equilK = equilibrium constant

% last modified on 4/5/18

    % gas constant
    idealR = 8.314E-3; % kJ/mol/K
    
    % species thermo
    [~, HsensCO, absSCO, dHformCO, ~] = thermoForCO(temperature);
    [~, HsensCH4, absSCH4, dHformCH4, ~] = thermoForCH4(temperature);
    [~, HsensH2O, absSH2O, dHformH2O, ~] = thermoForH2Ov(temperature);
    [~, HsensH2, absSH2, dHformH2, ~] = thermoForH2(temperature);
    
    % Enthalpy, entropy and free energy of reaction at 298K
    dHrxn298 = dHformCO + 3*dHformH2 - dHformCH4 - dHformH2O;
    
    % Enthalpy and free energy of reaction at specified T
    dHrxnT = -HsensCH4 - HsensH2O + dHrxn298 + HsensCO + 3*HsensH2;
    dSrxnT = absSCO + 3*absSH2 - absSCH4 - absSH2O;
    dGrxnT = dHrxnT - temperature*dSrxnT/1000;
    
    % Equilibrium constant at specified T
    equilK = exp(-dGrxnT/idealR/temperature);
end

