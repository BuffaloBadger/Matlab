function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForFeO(temperature)
%thermoForFeO Calculate thermochemical values for FeO
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 3/13/18
    enthalpyOfForm = -272.04;
    entropyAt298 = 60.75;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for FeO below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1650
        A = 45.7512;
        B = 18.78553;
        C = -5.952201;
        D = 0.852779;
        E = -0.081265;
        F = -286.7429;
        G = 110.3120;
        H = -272.0441;
    else
        display(' ')
        display('Thermo data for FeO(s) above 1650 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForFeO

