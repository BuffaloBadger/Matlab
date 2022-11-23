function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForNaCl(temperature)
%thermoForNaCl Calculate thermochemical values for NaCl
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 3/13/18
    enthalpyOfForm = -411.12;
    entropyAt298 = 72.11;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for NaCl below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1073
        A = 50.72389;
        B = 6.672267;
        C = -2.517167;
        D = 10.15934;
        E = -0.200675;
        F = -427.2115;
        G = 130.3973;
        H = -411.1203;
    else
        display(' ')
        display('Thermo data for NaCl above 1073 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForNaCl

