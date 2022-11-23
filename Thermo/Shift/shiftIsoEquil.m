function [equilMoleFract] = shiftIsoEquil(temperature, initMoleFract)
%shiftIsoEquil Calculate equilibrium composition for isothermal reaction
%   temperature = temperature in K
%   initMoleFract = vector of initial mole fractions where the elements are
%       1 - CO
%       2 - H2O
%       3 - CO2
%       4 - H2
%       5+ - any other species present in the system
%   equilMoleFract = vector of equilibrium mole fractions with the elements
%       in the same order as in initMoleFract

% last modified on 4/6/18
    
    % Basis of 1 mole total
    nInitCO = initMoleFract(1);
    nInitH2O = initMoleFract(2);
    nInitCO2 = initMoleFract(3);
    nInitH2 = initMoleFract(4);
    
    % guess the moles of the reagents at equilibrium
    nReactedGuess = 0.1*nInitCO;
    guess = [
        nInitCO - nReactedGuess;
        nInitH2O - nReactedGuess;
        nInitCO2 + nReactedGuess;
        nInitH2 + nReactedGuess;
    ];
    
    % calculate the moles of the reagents at equilibrium
    [equilMoles, flag, message] = solveATEs(@equilEqns, guess);
    
    % check for errors
    if flag < 0
        display(' ')
        display('Error solving the equilibrium equations:')
        display(message)
    end
    for i = 1:4
        if equilMoles(i) < 0
            display('Error: negative moles of a reagent')
        end
    end
    
    % return the equilibrium mole fractions of all species
    equilMoleFract = initMoleFract;
    equilMoleFract(1:4) = equilMoles(1:4);
    
    % internal function to calculate the residuals
    function residual = equilEqns(curMoleFract)
        % extract molar amounts as scalars
        nCO = curMoleFract(1);
        nH2O = curMoleFract(2);
        nCO2 = curMoleFract(3);
        nH2 = curMoleFract(4);
        
        % calculate the equilibrium constant
        equilK = shiftEquilK(temperature);
        
        % calculate the residuals
        residual = [equilK*nCO*nH2O - nH2*nCO2; % equilibrium eqn
            nInitCO + nInitCO2 - nCO - nCO2; % C balance
            2*nInitH2O + 2*nInitH2 - 2*nH2O - 2*nH2; % H balance
            nInitCO + nInitH2O + 2*nInitCO2 - nCO - nH2O - 2*nCO2 % O bal
        ];
    end
end

