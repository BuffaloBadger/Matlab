function [avg, hdi95] = analyzeMarkovChain(chain,plotTraces,plotHist)
%analyzeMarkovChain Analyze a Markov chain
%   For each chain the mean and 95% highest density interval are estimated.
%   Traces and histograms are optionally plotted.
%
%   last revised 11/13/17
%
%   avg - column vector of mean values of the parameters calculated from
%       the posterior pdf histograms
%   hdi95 - two-column matrix with lower limits of 95% hdi in the first
%       column and upper limits of 95% hdi in the second column again
%       calculated from the posterior pdf histograms
%   chain - matrix with Markov chains of each parameter as columns
	
    nVars = size(chain,2); % number of parameters to be analyzed
    % plot traces if requested
    if plotTraces
        for i = 1:nVars
            figure
            plot(chain(:,i))
            title(['Trace for Variable ',num2str(i,2)],'FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel('Sample Number', 'FontSize', 14)
            ylabel('Sample Value', 'FontSize', 14)
        end
    end % of plotTraces
    
    % calculate the means
    avg = mean(chain);
    
    hdi95 = [];
    % generate histogram data and plot the histograms, if requested
    for i = 1:nVars
        if plotHist
            % plot the posterior pdf histogram
            figure
            h = histogram(chain(:,i),'Normalization','pdf');
            % calculate the 95% HDI
            [min95hdi, max95hdi] = calcHDI(h.Values,h.BinEdges);
            hdi95 = [hdi95; min95hdi max95hdi];
            % add the mean and 95% hdi to the posterior pdf histogram
            hold on
            plot([avg(i) avg(i)],[0 1.1*max(h.Values)],'r','LineWidth',3)
            hold on
            plot([hdi95(i,1) hdi95(i,1)],[0 1.1*max(h.Values)],...
                'g','LineWidth',3)
            hold on
            plot([hdi95(i,2) hdi95(i,2)],[0 1.1*max(h.Values)],...
                'g','LineWidth',3)
            title(['Probability Density of Random Variable ',num2str(i,2)])
            set(gca, 'FontSize', 14);
            xlabel('Value of Random Variable', 'FontSize', 14)
            ylabel('Probability Density', 'FontSize', 14)
        else
            % get the histgram data
            [val, edges] = histcounts(chain(:,i),'Normalization','pdf');
            [min95hdi, max95hdi] = calcHDI(val,edges);
            hdi95 = [hdi95; min95hdi max95hdi];
        end
    end
    
    % Internal function to calculate the 95% hdi
    function [minHDI, maxHDI] = calcHDI(val,edges)
        binWidth = edges(2) - edges(1);
        % sort the histogram bars from largest to smallest
        [sortedVal, histIndex] = sort(val,'descend');
            
        % add histogram bars to the 95% hdi until the area equals 0.95
        area = 0.0;
        currentIndex = 1;
        while area < 0.95
            area = area + binWidth*sortedVal(currentIndex);
            currentIndex = currentIndex + 1;
        end
        % find the smallest and largest parameter values of the added
        % histogram bars
        minHDI = edges(histIndex(1));
        maxHDI = edges(histIndex(1) + 1);
        for j = 2:currentIndex
            if edges(histIndex(j)) < minHDI
                minHDI = edges(histIndex(j));
            end
            if edges(histIndex(j) + 1) > maxHDI
                maxHDI = edges(histIndex(j) + 1);
            end
        end
    end % of calcHDI
end % of analyzeMarkovChain

