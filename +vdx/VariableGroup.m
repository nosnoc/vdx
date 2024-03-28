classdef VariableGroup < handle &...
        matlab.mixin.indexing.RedefinesParen
    properties
        members
        
        %
        vector

        %
        solver
    end

    properties (Dependent)
    end
    
    methods(Access=public)
        function obj = VariableGroup(members, vector, solver)
            obj.members = members; % TODO(@anton) maybe make sure these are ok, i.e. all dims up to depth match and throw a good error
            obj.vector = vector;
            obj.solver = obj.solver;
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
                adj_ind = index_adjustment(index_op.Indices);
                symbolics = [];
                for v=members
                    var = obj.vector.(v);
                    depth = var.depth;
                    vsymbolics = cellfun(@(x) obj.vector.w(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                    symbolics = [symbolics;vsymbolics{:}];
                end
                out = squeeze(symbolics);
            else
                if index_op(2).Type == 'Dot'
                    switch(index_op(2).Name)
                      case "lb"
                        lb = [];
                        for v=members
                            var = obj.vector.(v);
                            depth = var.depth;
                            vlb = cellfun(@(x) obj.vector.lb(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            lb = [lb;vlb];
                        end
                        out = squeeze(lb);
                      case "ub"
                        ub = [];
                        for v=members
                            var = obj.vector.(v);
                            depth = var.depth;
                            vub = cellfun(@(x) obj.vector.ub(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            ub = [ub;vub];
                        end
                        out = squeeze(ub);
                      case "init"
                        init = [];
                        for v=members
                            var = obj.vector.(v);
                            depth = var.depth;
                            vinit = cellfun(@(x) obj.vector.init(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            init = [init;vinit];
                        end
                        out = squeeze(init);
                      case "res"
                        res = [];
                        for v=members
                            var = obj.vector.(v);
                            depth = var.depth;
                            vres = cellfun(@(x) obj.vector.res(x), var.indices(adj_ind{1:depth},1), 'uni', false);
                            res = [res;vres];
                        end
                        out = squeeze(res);
                      case "mult"
                        mult = [];
                        for v=members
                            var = obj.vector.(v);
                            depth = var.depth;
                            vmult = cellfun(@(x) obj.vector.mult(x), var.indices(adj_ind{1:depth},1), 'uni', false);
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

        function obj = parenAssign(obj,index_op,varargin)
            error('Assigning to a variable group is not supported')
        end

        function n = parenListLength(obj,index_op,ctx)
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
