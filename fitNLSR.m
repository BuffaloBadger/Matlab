function [beta, betaCI, rSquared] = fitNLSR(betaGuess, x, yMeas, calcY,...
    useRelErr)
%fitSR fit a model to single response data
%   1/9/23
%
%   betaGuess = column vector containing guesses for the parameter values
%   x = set variable matrix; x(i,j) is the value of set variable j in
%       experiment i
%   yMeas = measured response vector
%   calcY = function handle for a function that calculates the model-
%       predicted response vector given the parameter vector and the set
%       variable matrix
%   useRelErr = flag to use relative errors if true
%   beta = column vector of parameter values
%   betaCI = matrix containing the lower end of the 95% confidence interval
%       for beta in the first column and the upper end in the other column
%   rSquared = coefficient of determination for the fitted model
%
    if useRelErr
        weight = 1./(yMeas.^2);
        mdl = fitnlm(x,yMeas,calcY,betaGuess,'Weights',weight);
    else
        mdl = fitnlm(x,yMeas,calcY,betaGuess);
    end
    
    beta = mdl.Coefficients.Estimate;
    betaCI = coefCI(mdl);
    rSquared = mdl.Rsquared.Ordinary; % .Adjusted also available
end % of fitNLSR

