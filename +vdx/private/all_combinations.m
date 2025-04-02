function X = all_combinations(varargin)

    numSets = length(varargin);
    for i=1:numSets,
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
    for i=1:size(X,1)
        ixVect = cell(length(sizeThisSet),1);
        sz = flip(sizeThisSet);
        if length(sz) == 1
            sz = [sz,1];
        end
        [ixVect{:}] = ind2sub(sz,i);
        ixVect = flip([ixVect{:}]);
        vect = zeros(1, numSets);
        for jj=1:numSets
            vect(jj) = varargin{jj}(ixVect(jj));
        end
        X(i,:) = vect;
    end
end
