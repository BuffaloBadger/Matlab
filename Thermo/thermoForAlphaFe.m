function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForAlphaFe(temperature)
%thermoForAlphaFe Calculate thermochemical values for Fe
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 3/13/18
    enthalpyOfForm = 0.0;
    entropyAt298 = 27.31;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for alpha-Fe below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 700
        A = 18.42868;
        B = 24.64301;
        C = -8.91372;
        D = 9.664706;
        E = -0.012643;
        F = -6.573022;
        G = 42.51488;
        H = 0.0;
    elseif temperature < 1042
        A = -57767.65;
        B = 137919.7;
        C = -122773.2;
        D = 38682.42;
        E = 3993.08;
        F = 24078.67;
        G = -87364.01;
        H = 0.0;
    elseif temperature < 1100.
        A = -325.8859;
        B = 28.92876;
        C = 0;
        D = 0;
        E = 411.9629;
        F = 745.8231;
        G = 241.8766;
        H = 0.0;
    elseif temperature < 1809.
        A = -776.7387;
        B = 919.4005;
        C = -383.7184;
        D = 57.08148;
        E = 242.1369;
        F = 697.6234;
        G = -558.3674;
        H = 0.0;
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

