function [chain, pctAccepted] = sampleMarkovChain(target,initVal,...
    propFract,chainLen,burnFract,keepEveryOneIn,plotBurnIn)
%sampleMarkovChain generate a Markov chain
%   Metropolis-Hastings is used on one sampled variable at a time. The 
%   proposal pdf is a Normal distribution with the mean equal to the last 
%   value of the sampled variable and the standard deviation a constant
%   fraction of the mean. A fraction of the chain is removed from the 
%   beginning to allow burn-in, and the chain is thinned by using only
%   every nth sample.
%
%   last revised 12/5/17
%
%   chain - matrix with Markov chains of each parameter as columns
%   pctAccepted - percentage of the proposal values that were accepted
%
%   target - handle for a function that computes the probability density of
%       a sample of the quantities, to within a proportionality constant
%   initVal - column vector with initial values of quantities to be sampled
%   propFract - fraction of the proposal pdf mean to use as the proposal
%       pdf standard deviation
%   chainLen - length of the Markov chains prior to shortening for burn-in
%       and thinning
%   burnFract - fraction of the initial chains to be removed to allow
%       burn-in
%   keepEvery - chain elements with an index that is a multiple of this 
%       integer are retained during chain thinning
%   plotBurnIn - the portion of the chain removed as burn-in will be 
%       plotted if true

    % initialize
    nVars = length(initVal);
    chain = initVal.'; % place the initial values in the first row
    vcur = chain; % current set of values
    fcur = target(initVal);
    acceptCount = 0;

    % create the Markov chain
    for i = 2:chainLen
        for j = 1:nVars
            vprop = norminv(rand,vcur(j),propFract*vcur(j));
            while vprop < 0.0
                vprop = norminv(rand,vcur(j),propFract*vcur(j));
            end
            vcur(j) = vprop;
            fprop = target(vcur);
            a = fprop*normpdf(chain(i-1,j),vprop,propFract*vprop)/...
                fcur/normpdf(vprop,chain(i-1,j),propFract*chain(i-1,j));
            if a > rand
                fcur = fprop;
                acceptCount = acceptCount + 1;
            else
                vcur(j) = chain(i-1,j);
            end
        end
        chain = [chain; vcur];
    end
    
    % calculate the acceptance rate
    pctAccepted = 100.0*acceptCount/(nVars*(chainLen - 1));
	    
    % plot (if requested) and remove the burn-in
    newStart = fix(chainLen*burnFract) + 1;
    if plotBurnIn
        for i = 1:nVars
            figure
            plot(chain(1:newStart-1,i))
            title(['Burn-In Trace for Variable ',num2str(i,2)],...
                'FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel('Sample Number', 'FontSize', 14)
            ylabel('Sample Value', 'FontSize', 14)
        end
    end
    chain = chain(newStart:end,:);
	    
    % thin the chain
    chain = chain(1:keepEveryOneIn:end,:);
    
end % of sampleMarkovChain

