function [x, y] = solveNSBVODEs(xB1, xB2, yAvg, calcDerivs, calcResids)
%solveNSBVODEs solve a set of non-singular boundary value ODEs
%   8/3/21
%
%   xB1 = value of the independent variable at boundary 1
%   xB2 = value of the independent variable at boundary 2
%   yAvg = vector of average values of the dependent variables
%   calcDerivs = handle to a function that calculates the derivatives,
%       dy/dx given the a value of x values and a corresponding vector of
%       y values
%   calcResids = handle to a function that calculates the residuals at the
%       boundaries given a vector of y values at boundary 1 and a vector of
%       y values at boundary 2
%   x = vector of values of the independent variable at points starting at
%       boundary 1 and ending at boundary 2
%   y = matrix of values of the dependent variable j at point x(i)
%   dydx = vector of values of dy/dx at the points in x
%

	% Set up the initial mesh
	nMeshPoints = 100;
	x = linspace(xB1, xB2, nMeshPoints);

	% Create a structure containing the mesh and guesses
	solinit=bvpinit(x,yAvg);

	% Solve the odes
	result = bvp4c(calcDerivs,calcResids,solinit);
    
    x = result.x;
	y = result.y;

end % of solveNSBVODEs
