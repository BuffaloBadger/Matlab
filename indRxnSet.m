function [nInd,indepRxns,nuInd] = indRxnSet(nuFull)
%indRxnSet identifies a complete set of independent chemical reactions
%   last revised 11/8/17
%
%   nuFull(i,j) - stoichiometric coefficients for all species, i, and
%       reactions, j
%   nInd - the number of mathematically independent reactions
%   nuInd(i,j) - stoichiometric coefficients of the species, i, in the
%       mathematically independent reactions, j
%   indepRxns(j) - row vector containing the indices of the reactions in
%       the full stoichiometric coefficient matrix that are mathematically
%       independent
    nuFullTranspose = nuFull.';
    nInd = rank(nuFullTranspose);
    nTotal = size(nuFullTranspose,1);
    % the first reaction must be independent
    indepRxns = 1;
    indSet = nuFullTranspose(1,:);
    nAdded = 1;
    % check the rest
    for i = 2:nTotal
        temp = [C(i,:); indSet];
        if rank(temp) > nAdded
            indSet = temp;
            nAdded = nAdded + 1;
            indepRxns = [indepRxns, i];
        end
    end
    nuInd = indSet.';
end % of indRxnSet

