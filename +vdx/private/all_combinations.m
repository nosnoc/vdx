function X = all_combinations(varargin)

    numSets = length(varargin);
    for i = 1:numSets,
        thisSet = sort(varargin{i});
        if ~isequal(prod(size(thisSet)),length(thisSet)),
            error('All inputs must be vectors.')
        end
        if ~isnumeric(thisSet),
            error('All inputs must be numeric.')
        end
        sizeThisSet(i) = length(thisSet);
        varargin{i} = thisSet;
    end
    X = zeros(prod(sizeThisSet),numSets);
    for i = 1:size(X,1)
        
        ixVect = cell(length(sizeThisSet),1);
        [ixVect{:}] = ind2sub(sizeThisSet,i);
        ixVect = [ixVect{:}];
        for j = 1:numSets,
            X(i,j) = varargin{j}(ixVect(j));
        end
    end
end
