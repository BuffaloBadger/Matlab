function [ind, dep, flag] = solveIVODEs(ind0, dep0, fVar, fVal, ODEFcn, ...
    odesAreStiff, useMatrixForm, ODEMatrixFcn)
%solveIVODEs wrapper for solving initial value ODEs using ode45/ode15s
%   revised 11/10/17
%
%   ind0 = initial value of the independent variable
%   dep0 = column vector of initial values of the dependent variables
%   fVal = final value of either the independent variable or one of the 
%       dependent variables
%   fVar = integer equal to the index of the dependent variable whose final
%       value is equal to fVal, or equal to zero if fVal is the final value
%       of the independent variable
%   ODEFcn = handle to a function that returns a column vector containing 
%       the values of the derivatives at values of the independent and 
%       dependent variables that are passed to it as arguments
%   odesAreStiff = boolean that is true if the ODEs are stiff
%   useMatrixForm = boolean that is true if the ODEs are written as a
%       matrix equation: M*d(dep)/d(ind) = ODEFcn
%   ODEMatrixFcn = handle to a function that computes the matrix M if the
%       ODEs are written as a matrix equation: M*d(dep)/d(ind) = ODEFcn
%       or 0 if useMatrixForm is false. ODEMatrixFcn is passed the current
%       values of ind and dep and returns the matrix M.
%   ind = column vector containing values of the independent variable
%       spanning the range from the initial value to its final value
%   dep = matrix where each column contains the values of one of the
%       dependent variables at corresponding to the independent variable
%       values in ind
%   flag = integer flag where
%           1 indicates a solution was found
%          -1 indicates that the final value specified for the dependent
%               variable was not reached
%          -2 indicates that the step size may have been too large causing
%               the solution to be inaccurate
    if fVar == 0 % fVal is the final value of the independent variable
        flag = 1;
        if odesAreStiff
            if useMatrixForm
                options = odeset('Mass',ODEMatrixFcn);
                [ind,dep] = ode15s(ODEFcn,[ind0 fVal],dep0,options);
            else
                [ind,dep] = ode15s(ODEFcn,[ind0 fVal],dep0);
            end
        else % ODEs are not stiff
            if useMatrixForm
                options = odeset('Mass',ODEMatrixFcn);
                [ind,dep] = ode45(ODEFcn,[ind0 fVal],dep0,options);
            else
                [ind,dep] = ode45(ODEFcn,[ind0 fVal],dep0);
            end
        end
    else % fVal is the final value of dependent variable dep(fVar)
        options = odeset('Events',@stop);
        count = 0;
        flag = -1;
        indf = ind0 + 1.0;
        while (count < 10 && flag < 0)
            count = count + 1;
            % solve the odes
            if odesAreStiff
                if useMatrixForm
                    options = odeset(options,'Mass',ODEMatrixFcn);
                end
                [ind,dep] = ode15s(ODEFcn,[ind0 indf],dep0,options);
            else % ODEs are not stiff
                if useMatrixForm
                    options = odeset(options,'Mass',ODEMatrixFcn);
                end
                [ind,dep] = ode45(ODEFcn,[ind0 indf],dep0,options);
            end
            if ind(end) == indf % indf wasn't large enough
                indf = (ind0 + 1.0)*10^(count);
            elseif ind(end) < 0.1*indf % indf was too large
                indf = (ind0 + 1.0)/10^(count);
                flag = -2;
            else % indf was appropriate and the equations are solved
                flag = 1.0;
            end
        end
    end
    
    % Function that provides the integration stopping criterion
    function [stop_if_this_is_zero, isterminal, direction] = stop(...
            ind_cur,dep_cur)
        isterminal = 1;
        direction = 0;
        stop_if_this_is_zero = dep_cur(fVar) - fVal;
    end % of internal function stop

end % of solveIVODEs

