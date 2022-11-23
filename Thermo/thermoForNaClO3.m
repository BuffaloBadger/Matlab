function [heatCap, enthalpy, entropy, enthalpyOfForm, entropyAt298] =...
    thermoForNaClO3(temperature)
%thermoForNaClO3 Calculate thermochemical values for NaClO3
%   temperature = temperature in K
%   heatCap = constant pressure heat capacity in J/mol/K
%   enthalpy = standard enthalpy relative to enthalpy at 298.15 K in kJ/mol
%   entropy = standard entropy in J/mol/K
%   enthalpyOfForm = standard heat of formation at 298.15 K in kJ/mol
%   entropyAt298 = absolute standard entropy at 298.15 K in J/mol/K

% last modified on 3/13/18
    enthalpyOfForm = -365.4; % Wikipedia
    entropyAt298 = 129.7; % Wikipedia
    if temperature < 100
        display(' ')
        display('Thermo data for NaClO3 below 100 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    elseif temperature < 298
        % Shomate parameters calculated from Fransson's heat capacity
        % expression
        MW = 106.44;
        A = 0.0512*MW;
        B = 6.09*MW;
        C = -16.45*MW;
        D = 20.23*MW;
        E = 0;
        F = enthalpyOfForm - 0.298*A - B*(0.298)^2/2 - C*(0.298)^3/3 -...
            D*(0.298)^4/4;
        G = entropyAt298 - A*log(0.298) - 0.298*B - C*(0.298)^2/2 -...
            D*(0.298)^3/3;
        H = enthalpyOfForm;
    elseif temperature < 528
        % Shomate parameters calculated from Kelley's enthalpy expression
        % for the solid
        A = 4.184*13.07;
        B = 2*4.184*18.5;
        C = 0;
        D = 0;
        E = 0;
        F = -5.541*4.184 + enthalpyOfForm;
        G = entropyAt298 - A*log(0.298) - 0.298*B;
        H = enthalpyOfForm;
        %{
        % Shomate parameters calculated from Campbell's enthalpy expression
        % for the solid
        A = 4.184*10.92;
        B = 2*4.184*22.0;
        C = 0;
        D = 0;
        E = 0;
        F = -5.266*4.184 + enthalpyOfForm;
        G = entropyAt298 - A*log(0.298);
        H = enthalpyOfForm;
        %}
    elseif temperature < 600
        % Shomate parameters calculated from Kelley's enthalpy expression
        % for the liquid
        A = 4.184*31.8;
        B = 0;
        C = 0;
        D = 0;
        E = 0;
        F = -4.87*4.184 + enthalpyOfForm;
        G = entropyAt298 - A*log(0.298);
        H = enthalpyOfForm;
        %{
        % Shomate parameters calculated from Campbell's enthalpy expression
        % for the liquid
        A = 4.184*32.05;
        B = 0;
        C = 0;
        D = 0;
        E = 0;
        F = -5.171*4.184 + enthalpyOfForm;
        G = entropyAt298 - A*log(0.298);
        H = enthalpyOfForm;
        %}
    else
        display(' ')
        display('Thermo data for NaClO3 above 600 K are not available')
        heatCap = NaN;
        enthalpy = NaN;
        entropy = NaN;
        return
    end
    t = temperature/1000.;
    heatCap = A + B*t + C*t^2 + D*t^3 + E/t^2;
    enthalpy = A*t + B/2*t^2 + C/3*t^3 + D/4*t^4 - E/t + F - H;
    entropy = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/2/t^2 + G;
end % of thermoForNaClO3

