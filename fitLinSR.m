function [beta, betaCI, rSquared] = fitLinSR(x, y, addIntercept, ...
    interceptOnly, makeModelPlot, makeResidualsPlots)
%fitLinSR fit a linear model to single response data
%
%   x = independent variable matrix; rows are different experiments, 
%       columns are different independent variables (x1, x2, etc.)
%   y = measured response variable column vector
%   addIntercept = boolean flag to include an intercept if true
%   interceptOnly = boolean flag to indicate that the model y = x + b
%       should be used
%   makeModelPlot = boolean flag to generate a model or parity plot
%   makeResidualsPlots = boolean flag to generate residuals plots
%   beta = parameter column vector
%       if an intercept was added, it is the last element
%   betaCI = parameter 95% confidence interval column vector
%   rSquared = coefficient of determination for the fitted model
%
%   last revised 1/14/20
    
    % initialize
    nData = length(y);
    nX = size(x,2);
    
    if interceptOnly % the model is y = x + beta
        % estimate the intercept
        beta = sum(y - x)/nData;
        
        % calculate the model-predicted responses and residuals
        yModel = x + beta;
        eps = y - yModel;
        
        % calculate the 95% confidence interval
        tCrit = tinv(0.975,(nData - 1));
        betaCI = tCrit*sqrt(sum((y - yModel).^2)/(nData - 1));
    else % the model is y = x*beta
        if addIntercept % add a column containing 1's to the end of x
            x = [x, ones(nData,1)];
        end
        
        % calculate the model parameters
        beta = inv(x.'*x)*x.'*y;
        
        % calculate the model-predicted responses
        yModel = x*beta;
        
        % calculate the parameter 95% confidence intervals
        nParams = length(beta);
        eps = y - x*beta;
        varEps = 1/(nData - nParams)*(eps.'*eps);
        covBeta = varEps*inv(x.'*x);
        tCrit = tinv(0.975,(nData - nParams));
        betaCI = [];
        for i = 1:nParams
            ci = sqrt(covBeta(i,i))*tCrit;
            betaCI = [betaCI ; ci];
        end
    end % of calculating the model parameters and their confidence 
        % intervals
    
    % calculate the coefficient of determination
    ybar = sum(y)/nData;
    rSquared = (sum((y-ybar).^2) - sum((y-yModel).^2))/sum((y-ybar).^2);

    if makeModelPlot
        if nX == 1 % make a model plot
            figure;
            plot(x(:,1),yModel,'r',x(:,1),y,'ok','LineWidth',2)
            title('Model Plot','FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel('Set Variable x', 'FontSize', 14)
            ylabel('Response Variable y','FontSize', 14)
            legend({'Model Prediction','Measured Value'}, 'Location',...
                'northwest', 'FontSize',14)
        else % make a parity plot
            figure;
            plot(y,y,'r',y,yModel,'ok','LineWidth',2)
            title('Parity Plot','FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel('Measured Response', 'FontSize', 14)
            ylabel('Predicted Response', 'FontSize', 14)
            legend({'Perfect Model/Data','Actual Model/Data'},...
                'Location','northwest','FontSize',14)
        end
    end % of making the model/parity plot
    
    if makeResidualsPlots
        for n = 1:nX
            figure;
            plot(x(:,n),zeros(nData,1),'r',x(:,n),eps,'ok','LineWidth',2)
            title(['Residuals Plot vs. x(',num2str(n),')'],'FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel(['Set variable x(',num2str(n),')'], 'FontSize', 14)
            ylabel('Residual', 'FontSize', 14)
        end
    end % of making the parity plots
end % of fitLinSR

