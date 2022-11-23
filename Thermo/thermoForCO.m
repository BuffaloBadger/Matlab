function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForCO(temperature)
%thermoForCO Calculate thermochemical values for CO
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    % Shomate equation parameters from the NIST Chemistry Web Book
    enthalpyOfForm = -110.53;
    entropyAt298 = 197.66;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for CO below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1300
        A = 25.56759;
        B = 6.096130;
        C = 4.054656;
        D = -2.671301;
        E = 0.131021;
        F = -118.0089;
        G = 227.3665;
        H = -110.5271;
    elseif temperature < 6000
        A = 35.15070;
        B = 1.300095;
        C = -0.205921;
        D = 0.013550;
        E = -3.282780;
        F = -127.8375;
        G = 231.7120;
        H = -110.5271;
    else
        display(' ')
        display('Thermo data for CO above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForCO

