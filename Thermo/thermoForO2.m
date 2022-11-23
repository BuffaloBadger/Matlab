function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForO2(temperature)
%thermoForO2 Calculate thermochemical values for O2
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 5/23/18
    enthalpyOfForm = 0.0;
    entropyAt298 = 205.15;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 100
        disp(' ')
        disp('Thermo data for O2 below 100 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 700
        A = 31.32234;
        B = -20.23531;
        C = 57.86644;
        D = -36.50624;
        E = -0.007374;
        F = -8.903471;
        G = 246.7945;
        H = 0.0;
    elseif temperature < 2000
        A = 30.03235;
        B = 8.772972;
        C = -3.988133;
        D = 0.788313;
        E = -0.741599;
        F = -11.32468;
        G = 236.1663;
        H = 0.0;
    elseif temperature < 6000
        A = 20.91111;
        B = 10.72071;
        C = -2.020498;
        D = 0.146449;
        E = 9.245722;
        F = 5.337651;
        G = 237.6185;
        H = 0.0;
    else
        disp(' ')
        disp('Thermo data for O2 above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForO2

