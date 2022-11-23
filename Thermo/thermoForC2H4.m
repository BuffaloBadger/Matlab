function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForC2H4(temperature)
%thermoForC2H4 Calculate thermochemical values for C2H4
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/7/20
    enthalpyOfForm = 52.47;
    entropyAt298 = 219.32;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for C2H4 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1200
        A = -6.387880;
        B = 184.4019;
        C = -112.9718;
        D = 28.49593;
        E = 0.315540;
        F = 48.17332;
        G = 163.1568;
        H = 52.46694;
    elseif temperature < 6000
        A = 106.5104;
        B = 13.73260;
        C = -2.628481;
        D = 0.174595;
        E = -26.14469;
        F = -35.36237;
        G = 275.0424;
        H = 52.46694;
    else
        display(' ')
        display('Thermo data for C2H4 above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForCH4

