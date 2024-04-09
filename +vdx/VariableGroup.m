classdef VariableGroup < handle &...
        matlab.mixin.indexing.RedefinesParen & ...
        matlab.mixin.indexing.RedefinesBrace
    properties
        members

        % Index rules must each be an object that takes an index (or many indices as cells) and a vdx.Variable and return an n
        % this should be a dictionary.
        index_rules
    end

    properties (Dependent)
    end
    
    methods(Access=public)
        function obj = VariableGroup(members, varargin)
            p = inputParser;
            addOptional(p, 'index_rules', {});
            parse(p, varargin{:});
            
            obj.members = members; % TODO(@anton) maybe make sure these are ok, i.e. all dims up to depth match and throw a good error
            obj.index_rules = p.Results.index_rules;
            
            % Populate index rules with identity.
            for ii=1:length(obj.members)
                if ii > length(obj.index_rules) || isempty(obj.index_rules{ii})
                    obj.index_rules{ii} = @vdx.indexing.identity;
                end
            end
        end
        
        function out = cat(dim,varargin)
            error('Concatenation not supported')
        end

        function varargout = size(obj,varargin)
            varargout = {size(obj.indices, varargin{:})};
            %TODO(anton) needs to return correct values for varargin
        end
        
        function ind = end(obj,k,n)
            sz = size(obj.indices);
            ind = sz(k)-1;
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            members = obj.members;
            if ~is_index_scalar(index_op(1).Indices)
                error('VariableGroups only permit scalar indexing for now.')
            end
            if isscalar(index_op)
                symbolics = [];
                for ii=1:length(members)
                    var = members{ii};
                    depth = var.depth;
                    index_rule = obj.index_rules{ii};
                    adj_ind = index_rule(index_op.Indices, var);
                    vsymbolics = cellfun(@(x) var.vector.sym(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                    symbolics = [symbolics;vsymbolics{:}];
                end
                out = squeeze(symbolics);
            else
                if index_op(2).Type == 'Dot'
                    switch(index_op(2).Name)
                      case "lb"
                        lb = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vlb = cellfun(@(x) var.vector.lb(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            lb = [lb;vlb];
                        end
                        out = squeeze(lb);
                      case "ub"
                        ub = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vub = cellfun(@(x) var.vector.ub(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            ub = [ub;vub];
                        end
                        out = squeeze(ub);
                      case "init"
                        init = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vinit = cellfun(@(x) var.vector.init(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            init = [init;vinit];
                        end
                        out = squeeze(init);
                      case "res"
                        res = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vres = cellfun(@(x) var.vector.res(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            res = [res;vres];
                        end
                        out = squeeze(res);
                      case "mult"
                        mult = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vmult = cellfun(@(x) var.vector.mult(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            mult = [mult;vmult];
                        end
                        out = squeeze(mult);
                      otherwise
                        error('vdx only supports getting lb, ub, init, res, or mult for a variable group via dot indexing');
                    end
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here.
                end
            end
            if isscalar(out)
                varargout{1} = out{1};
            else
                out = permute(out, ndims(out):-1:1);
                varargout{1} = [out{:}];
            end
        end

        function varargout = braceReference(obj, index_op)
            members = obj.members;
            if ~is_index_scalar(index_op(1).Indices)
                error('VariableGroups only permit scalar indexing for now.')
            end
            if isscalar(index_op)
                symbolics = [];
                for ii=1:length(members)
                    var = members{ii};
                    depth = var.depth;
                    index_rule = obj.index_rules{ii};
                    adj_ind = index_rule(index_op.Indices(1:depth), var);
                    vsymbolics = cellfun(@(x) var.vector.sym(x), var.indices(adj_ind{:},1), 'uni', false);
                    symbolics = [symbolics;vsymbolics];
                end
                out = squeeze(symbolics);
            else
                if index_op(2).Type == 'Dot'
                    switch(index_op(2).Name)
                      case "lb"
                        lb = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vlb = cellfun(@(x) var.vector.lb(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            lb = [lb;vlb];
                        end
                        out = squeeze(lb);
                      case "ub"
                        ub = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vub = cellfun(@(x) var.vector.ub(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            ub = [ub;vub];
                        end
                        out = squeeze(ub);
                      case "init"
                        init = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vinit = cellfun(@(x) var.vector.init(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            init = [init;vinit];
                        end
                        out = squeeze(init);
                      case "res"
                        res = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vres = cellfun(@(x) var.vector.res(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            res = [res;vres];
                        end
                        out = squeeze(res);
                      case "mult"
                        mult = [];
                        for ii=1:length(members)
                            var = members{ii};
                            depth = var.depth;
                            index_rule = obj.index_rules{ii};
                            adj_ind = index_rule(index_op.Indices, var);
                            vmult = cellfun(@(x) var.vector.mult(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            mult = [mult;vmult];
                        end
                        out = squeeze(mult);
                      otherwise
                        error('vdx only supports getting lb, ub, init, res, or mult for a variable group via dot indexing');
                    end
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here.
                end
            end
            varargout{1} = out;
        end

        function obj = parenAssign(obj,index_op,varargin)
            error('Assigning to a variable group is not supported')
        end

        function n = parenListLength(obj,index_op,ctx)
           n = 1;
        end

        function obj = braceAssign(obj,index_op,varargin)
            error('Assigning to a variable group is not supported')
        end

        function n = braceListLength(obj,index_op,ctx)
           n = 1;
        end

        function obj = parenDelete(obj,index_op)
            error('Deletion of symbolics is not supported through the variable view')
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            obj = VariableGroup([],[],[]);
        end
    end
end

function res = is_index_scalar(index)
    res = all(cellfun(@(x) isscalar(x) & ~ischar(x), index));
end
