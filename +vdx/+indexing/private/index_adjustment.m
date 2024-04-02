function indices = index_adjustment(indices)
% INDEX_ADJUSTMENT  Adjusts indices to match 0 indexing convention.

    for ii=length(indices):-1:1
        % we go backwards because in the future we may want to also handle negative indices as affecting higher level index values
        if ~ischar(indices{ii})
            indices{ii} = indices{ii} + 1;
        end
    end
end
