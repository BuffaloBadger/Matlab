function sample = sampleMonteCarlo(probDensFunc, lowerLimit, guess)
%sampleMonteCarlo draw a sample from distribution
%   last revised 12/4/17
%
%   probDensFunc = handle to a function that takes a vector of values as an 
%       argument and returns a vector of the corresponding the probability 
%       densities of those values
%   lowerLimit = the smallest possible value of the distributed quantity
%   guess = a guess for the value (used in the numerical of the value) 
%
    u = rand;
    [sample, flag, message] = solveATEs(@cumDensFunc, guess);
    if flag < 0
        display(' ')
        display(['Error sampling the distribution: ',message])
    end
    
    % Function that evaluates the residual of the cumulative density
    function residual = cumDensFunc(x)
        residual = u - integral(probDensFunc, lowerLimit, x);
    end % of cumDensFunc
end % of sampleMonteCarlo

