function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForH2Ol(temperature)
%thermoForH2Ol Calculate thermochemical values for H2O liquid
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    enthalpyOfForm = -285.83;
    entropyAt298 = 69.95;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for H2O liquid below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 500
        A = -203.6060;
        B = 1523.290;
        C = -3196.413;
        D = 2474.455;
        E = 3.855326;
        F = -256.5478;
        G = -488.7163;
        H = -285.8304;
    else
        display(' ')
        display('Thermo data for H2O liquid above 500 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        enthalpyOfForm = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForH2Ol

