function [beta, betaCI, rSquared] = fitSR(betaGuess, x, yMeas, calcY,...
    useRelErr)
%fitSR fit a model to single response data
%   8/3/21
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
%   betaCI = column vector containing the 95% confidence intervals for the 
%       parameters in beta
%   rSquared = coefficient of determination for the fitted model
%
    
    % Calculate the maximum likelihood parameter estimates and their 95%
    % confidence intervals
    if useRelErr
        weight = 1./yMeas;
    else
        weight = ones(size(yMeas,1),1);
    end
    [beta,~,~,sigma,~] = nlinfit(x,yMeas,calcY,betaGuess,'Weights',weight);
    eps = (yMeas - calcY(beta,x));
    confint = nlparci(beta,eps,'covar',sigma);
    betaCI = (confint(:,2) - confint(:,1))/2.0;
    
    % Calculate the coefficient of determination
    rSquared = 1.0 - (sum(eps.^2))/(sum((yMeas-mean(yMeas)).^2));
    
end % of fitNonLinSR

