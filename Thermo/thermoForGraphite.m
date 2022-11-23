function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForGraphite(temperature)
%thermoForAlphaFe Calculate thermochemical values for Fe
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 8/6/19
    enthalpyOfForm = 0.0;
    entropyAt298 = 5.6;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 273
        display(' ')
        display('Thermo data for alpha-Fe below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 3650
        A = 8.12;
        B = 0.0;
        C = 0.0;
        D = 0.0;
        E = 0.0;
        F = 0.0;
        G = entropyAt298 - A*log(0.298);
        H = A*0.298;
    else
        display(' ')
        display('Thermo data for alpha Fe above 1809 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForAlphaFe

