function out_indices = next(indices, variable)
    n_idx = size(indices(1), 1);
    idx_depth = variable.depth;
    out_indices = indices;
    varsize = size(variable);
    max_ind = cumprod(varsize);
    adj_indices = index_adjustment(indices);

    inds = sub2ind(flip(varsize), adj_indices{:});

    for ii=1:length(inds)
        ind = inds(ii)+1;
        sub = cell(1, idx_depth);
        [sub{:}] = ind2sub(flip(varsize), ind);
        sub = flip(sub);
        
        while(isempty(variable.indices{sub{:}}))
            ind = ind+1;
            if(ind > max_ind)
                error("failed to find a next value")
            end
            [sub{:}] = ind2sub(flip(varsize), ind);
            sub = flip(sub);
            sub
            ind
            variable.indices{sub{:}}
        end
        inds(ii) = ind;
    end

    out_indices = cell(1, idx_depth);
    [out_indices{:}] = ind2sub(flip(varsize), inds);
    out_indices=flip(out_indices);
end
