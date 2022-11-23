function vaporPressure = vaporPressForH2O(temperature)
%vaporPressForH2O calculate the vapor pressure of water
%   temperature = temperature in K
%   vaporPressure = vapor pressure in atm
%
%   Antoine equation data from the NIST Chemistry Webbook

% last modified 4/20/18
    if temperature < 255.9
        display(' ')
        display(['Vapor pressure data are not available for H2O below',...
            '255.9 K'])
        vaporPressure = NaN;
        return
    elseif temperature < 373.0 % use the Stull parameters
        vaporPressure = 10^(0.98692*(4.6543 - (1435.264/(temperature -...
            64.848))));
    elseif temperature < 379.0 % use the average of the Stull and Liu
        % parameters
        vaporPressureStull = 10^(0.98692*(4.6543 - (1435.264/(temperature -...
            64.848))));
        vaporPressureLiu = 10^(0.98692*(3.55959 - (643.748/(temperature -...
            198.043))));
        vaporPressure = (vaporPressureStull + vaporPressureLiu)/2.0;
    elseif temperature < 573.0 % use the Liu parameters
        vaporPressure = 10^(0.98692*(3.55959 - (643.748/(temperature -...
            198.043))));
    else
        display(' ')
        display(['Vapor pressure data are not available for H2O above',...
            '573 K'])
        vaporPressure = NaN;
    end
end

