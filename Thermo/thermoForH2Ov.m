function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForH2Ov(temperature)
%thermoForH2Ov Calculate thermochemical values for H2O vapor
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 2/14/18
    enthalpyOfForm = -241.83;
    entropyAt298 = 188.84;
    if temperature < 500
        % A through D from Reid Prausnizt and Sherwood, converted from cal
        % to J and T to t
        A = 7.701*4.184;
        B = 0.4595*4.184;
        C = 2.521*4.184;
        D = -0.859*4.184;
        % E set to zero since Reid, Prausnitz and Sherwood don't use it
        E = 0.0;
        % F and G calculated
        F = enthalpyOfForm - A*0.298 - B/2*0.298^2 - C/3*0.298^3 -...
            D/4*0.298^4;
        G = entropyAt298 - A*log(0.298) - B*0.298 - C/2*0.298^2 -...
            D/3*0.298^3;
        H = enthalpyOfForm;
    % Shomate equation parameters from the NIST Chemistry Web Book
    elseif temperature < 1700
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
        F = -250.8810;
        G = 223.3967;
        H = -241.8264;
    elseif temperature < 6000
        A = 41.96426;
        B = 8.622053;
        C = -1.499780;
        D = 0.098119;
        E = -11.15764;
        F = -272.1797;
        G = 219.7809;
        H = -241.8264;
    else
        display(' ')
        display('Thermo data for H2O vapor above 6000 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForH2Ov

