function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForCH4(temperature)
%thermoForCH4 Calculate thermochemical values for CH4
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    enthalpyOfForm = -74.6;
    entropyAt298 = 186.25;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for CH4 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1300
        A = -0.703029;
        B = 108.4773;
        C = -42.52157;
        D = 5.862788;
        E = 0.678565;
        F = -76.84376;
        G = 158.7163;
        H = -74.8731;
    elseif temperature < 6000
        A = 85.81217;
        B = 11.26467;
        C = -2.114146;
        D = 0.138190;
        E = -26.42221;
        F = -153.5327;
        G = 224.4143;
        H = -74.8731;
    else
        display(' ')
        display('Thermo data for CH4 above 6000 K are not available')
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

