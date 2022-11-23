function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForCO2(temperature)
%thermoForCO2 Calculate thermochemical values for CO2
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    enthalpyOfForm = -393.52;
    entropyAt298 = 213.79;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Therom data for CO2 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1200
        A = 24.99735;
        B = 55.18696;
        C = -33.69137;
        D = 7.948387;
        E = -0.136638;
        F = -403.6075;
        G = 228.2431;
        H = -393.5224;
    elseif temperature < 6000
        A = 58.16639;
        B = 2.720074;
        C = -0.492289;
        D = 0.038844;
        E = -6.447293;
        F = -425.9186;
        G = 263.6125;
        H = -393.5224;
    else
        display(' ')
        display('Thermo data for CO2 above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForCO2

