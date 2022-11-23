function [dryMoleFract,totalMoleFract] = steamReformEquil(temperature,...
    pressure,steamToMethane)
%steamReformEquil Equilibrium composition from steam reforming
%   temperature = temperature in K
%   pressure = pressure in atm
%   steamToMethane = molar ratio of steam to CH4 in the feed
%   dryMoleFract = dry gas mole fractions
%       1 - H2
%       2 - CO
%       3 - CO2
%       4 - CH4
%   totalMoleFract = wet gas mole fractions
%       1 - H2
%       2 - CO
%       3 - CO2
%       4 - CH4
%       5 - H2O

% last modified on 4/5/18
    
    % Basis: 1 mol CH4
    n0CH4 = 1;
    n0H2O = steamToMethane;
    
    % guess the moles of the species at equilibrium
    guess = [3.0; % moles of H2
        0.07; % moles CO
        0.92; % moles CO2
        0.01; % moles CH4
        n0H2O - 1]; % moles H2O
    
    % calculate the moles of the species at equilibrium
    [n, flag, message] = solveATEs(@evaluateResiduals, guess);
    
    % check for errors
    if flag < 0
        display(' ')
        display('Error solving the equilibrium equations:')
        display(message)
    end
    for i = 1:5
        if n(i) < 0.0
            display(' ')
            display('Calculated equilibrium moles are negative')
        end
    end
    
    % calculate the total and dry mole fractions
    nTotal = sum(n);
    totalMoleFract = n/nTotal;
    nDry = n(1:4);
    nTotalDry = sum(nDry);
    dryMoleFract = nDry/nTotalDry;
    
    % internal function to evaluate the equilibrium expressions as
    % residuals
    function residual = evaluateResiduals(nCurrentGuess)
        % extract molar amounts as scalars
        nH2 = nCurrentGuess(1);
        nCO = nCurrentGuess(2);
        nCO2 = nCurrentGuess(3);
        nCH4 = nCurrentGuess(4);
        nH2O = nCurrentGuess(5);
        nTot = sum(nCurrentGuess);
        
        % calculate the equilibrium constants
        Kwgs = shiftEquilK(temperature);
        Ksr = steamReformingEquilK(temperature);
        
        % calculate the residuals
        residual = [Ksr*nCH4*nH2O*nTot^2 - nCO*nH2^3*pressure^2;
            Kwgs*nCO*nH2O - nCO2*nH2;
            n0CH4 - nCO - nCO2 - nCH4;
            n0H2O - nH2O - nCO - 2*nCO2;
            4*n0CH4 + 2*n0H2O - 4*nCH4 - 2*nH2O - 2*nH2
        ];
    end
end

