function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForH2(temperature)
%thermoForH2 Calculate thermochemical values for H2
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    enthalpyOfForm = 0.0;
    entropyAt298 = 130.68;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for H2 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1000
        A = 33.066178;
        B = -11.363417;
        C = 11.432816;
        D = -2.772874;
        E = -0.158558;
        F = -9.980797;
        G = 172.707974;
        H = 0.0;
    elseif temperature < 2500
        A = 18.563083;
        B = 12.257357;
        C = -2.859786;
        D = 0.268238;
        E = 1.977990;
        F = -1.147438;
        G = 156.288133;
        H = 0.0;
    elseif temperature < 6000
        A = 43.413560;
        B = -4.293079;
        C = 1.272428;
        D = -0.096876;
        E = -20.533862;
        F = -38.515158;
        G = 162.081354;
        H = 0.0;
    else
        display(' ')
        display('Thermo data for H2 above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForH2

