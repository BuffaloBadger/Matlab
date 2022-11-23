function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForC2H6(temperature)
%thermoForC2H6 Calculate thermochemical values for C2H6
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K


% last modified on 2/7/20
    enthalpyOfForm = -83.820; % Smith, Van Ness, Abbott & Swihart
    entropyAt298 = 229.6; % Wikipedia
    % Shomate equation parameters calculated from expressions in Smith, 
    % Van Ness, Abbott & Swihart
    if temperature < 298
        display(' ')
        display('Thermo data for C2H6 below 298 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 1500
        R = 8.3145;
        A = 1.131*R;
        B = 19.225*R;
        C = -5.561*R;
        D = 0.0;
        E = 0.0;
        F = enthalpyOfForm - 0.298*A - B*(0.298)^2/2 - C*(0.298)^3/3;
        G = entropyAt298 - A*log(0.298) - 0.298*B - C*(0.298)^2/2;
        H = enthalpyOfForm;
    else
        display(' ')
        display('Thermo data for C2H6 above 1500 K are not available')
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

