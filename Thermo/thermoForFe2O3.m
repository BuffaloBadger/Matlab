function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForFe2O3(temperature)
%thermoForFe2O3 Calculate thermochemical values for Fe2)3
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 3/13/18
    enthalpyOfForm = -825.50;
    entropyAt298 = 87.28;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for Fe2O3 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 950
        A = 93.43834;
        B = 108.3577;
        C = -50.86447;
        D = 25.58683;
        E = -1.61133;
        F = -863.2094;
        G = 161.0719;
        H = -825.5032;
    elseif temperature < 1050
        A = 150.624;
        B = 0;
        C = 0;
        D = 0;
        E = 0;
        F = -875.6066;
        G = 252.8814;
        H = -825.5032;
    elseif temperature < 2500
        A = 110.9362;
        B = 32.04714;
        C = -9.192333;
        D = 0.901506;
        E = 5.433677;
        F = -843.1471;
        G = 228.3548;
        H = -825.5032;
    else
        display(' ')
        display('Thermo data for Fe2O3 above 2500 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForFe2O3

