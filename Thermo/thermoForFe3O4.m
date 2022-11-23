function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForFe3O4(temperature)
%thermoForFe3O4 Calculate thermochemical values for Fe3O4
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

    % last modified 3/13/18
    enthalpyOfForm = -1120.89;
    entropyAt298 = 145.2;
    % Shomate equation parameters from the NIST Chemistry Web Book
    if temperature < 298
        display(' ')
        display('Thermo data for Fe3O4 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 900
        A = 104.2096;
        B = 178.5108;
        C = 10.6151;
        D = 1.132534;
        E = -0.994202;
        F = -1163.336;
        G = 212.0585;
        H = -1120.894;
    elseif temperature < 3000
        A = 200.832;
        B = 1.586435E-7;
        C = -6.661682E-8;
        D = 9.452452E-9;
        E = 3.186020E-8;
        F = -1174.135;
        G = 388.0790;
        H = -1120.894;
    else
        display(' ')
        display('Thermo data for Fe3O4 above 3000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForFe3O4

