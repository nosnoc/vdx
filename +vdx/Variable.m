classdef Variable < handle &...
        matlab.mixin.indexing.RedefinesBrace &...
        matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.Copyable
    properties
        %
        indices

        %
        vector

        %
        solver
    end

    methods(Access=public)
        function obj = Variable(vector, solver)
            obj.indices = cell(0,1);
            obj.vector = vector;
            obj.solver = obj.solver;
        end
        
        function out = cat(dim,varargin)
            error('Concatenation not supported')
        end

        function varargout = size(obj,varargin)
            varargout = size(obj.indices, varargin{:});
            %TODO(anton) needs to return correct values for varargin
        end
        
        function ind = end(obj,k,n)
            sz = size(obj.indices);
            ind = sz(k)-1;
        end

        function output = print(obj, varargin)
        % TODO(@anton)
            w = false;
            lb = false;
            ub = false;
            init = false;
            res = false;
            mult = false;
            if isempty(varargin)
                w = true;
                lb = true;
                ub = true;
                init = true;
                res = true;
                mult = true;
            else
                if any(ismember(lower(varargin), 'w'))
                    w = true;
                end
                if any(ismember(lower(varargin), 'lb'))
                    lb = true;
                end
                if any(ismember(lower(varargin), 'ub'))
                    ub = true;
                end
                if any(ismember(lower(varargin), 'init'))
                    init = true;
                end
                if any(ismember(lower(varargin), 'res'))
                    res = true;
                end
                if any(ismember(lower(varargin), 'mult'))
                    mult = true;
                end
            end

            % Generate header
            header = 'i\t\t';
            if lb
                header = [header 'lb\t\t'];
            end
            if ub
                header = [header 'ub\t\t'];
            end
            if init
                header = [header 'init\t\t'];
            end
            if res
                header = [header 'res\t\t'];
            end
            if mult
                header = [header 'mult\t\t'];
            end
            if w
                header = [header 'sym\t\t'];
            end
            header = [header '\n'];

            % iterate over all requested values
            indices = sort([obj.indices{:}]);
            output = header;
            for ii=indices
                pline = [num2str(ii) '\t\t'];
                if lb
                    pline = [pline sprintf('%-8.5g\t', obj.vector.lb(ii))];
                end
                if ub
                    pline = [pline sprintf('%-8.5g\t', obj.vector.ub(ii))];
                end
                if init
                    pline = [pline sprintf('%-8.5g\t', obj.vector.init(ii))];
                end
                if res
                    pline = [pline sprintf('%-8.5g\t', obj.vector.res(ii))];
                end
                if mult
                    pline = [pline sprintf('%-8.5g\t', obj.vector.mult(ii))];
                end
                if w
                    pline = [pline char(formattedDisplayText(obj.vector.w(ii)))];
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            fprintf(output);
        end
    end

    methods (Access=protected)
        function varargout = parenReference(obj, index_op)
            if isscalar(index_op)
                % TODO(@anton) Decide whether we squeeze, or concatenate with sorted indices.
                %              This is in my opinion purely a decision that should be made and stuck to.
                adj_ind = index_adjustment(index_op.Indices);
                symbolics = cellfun(@(x) obj.vector.w(x), obj.indices(adj_ind{:}), 'uni', false);
                out = squeeze(symbolics);
            else
                if index_op(2).Type == 'Dot'
                    switch(index_op(2).Name)
                      case "lb"
                        lb = cellfun(@(x) obj.vector.lb(x), obj.indices, 'uni', false);
                        adj_ind = index_adjustment(index_op(1).Indices);
                        out = squeeze(lb(adj_ind{:}));
                      case "ub"
                        ub = cellfun(@(x) obj.vector.ub(x), obj.indices, 'uni', false);
                        adj_ind = index_adjustment(index_op(1).Indices);
                        out = squeeze(ub(adj_ind{:}));
                      case "init"
                        init = cellfun(@(x) obj.vector.init(x), obj.indices, 'uni', false);
                        adj_ind = index_adjustment(index_op(1).Indices);
                        out = squeeze(init(adj_ind{:}));
                      case "res"
                        res = cellfun(@(x) obj.vector.res(x), obj.indices, 'uni', false);
                        adj_ind = index_adjustment(index_op(1).Indices);
                        out = squeeze(res(adj_ind{:}));
                      case "mult"
                        mult = cellfun(@(x) obj.vector.mult(x), obj.indices, 'uni', false);
                        adj_ind = index_adjustment(index_op(1).Indices);
                        out = squeeze(mult(adj_ind{:}));
                      otherwise
                        error('vdx only supports getting lb, ub, init, res, or mult for a variable via dot indexing');
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
            if isscalar(index_op)
                % get the cell array of args (x,x0,lbx,ubx)
                arg = varargin{1};
                indices = obj.vector.add_variable(arg{:});
                adj_ind = index_adjustment(index_op.Indices);
                obj.indices{adj_ind{:}} = indices;
            else
                if index_op(2).Type == 'Dot'
                    if all(cellfun(@(x) isscalar(x) & ~ischar(x), index_op(1).Indices))
                        switch(index_op(2).Name)
                          case "lb"
                            adj_ind = index_adjustment(index_op(1).Indices);
                            obj.vector.lb(obj.indices{adj_ind{:}}) = varargin{1};
                          case "ub"
                            adj_ind = index_adjustment(index_op(1).Indices);
                            obj.vector.ub(obj.indices{adj_ind{:}}) = varargin{1};
                          case "init"
                            adj_ind = index_adjustment(index_op(1).Indices);
                            obj.vector.init(obj.indices{adj_ind{:}}) = varargin{1};
                          otherwise
                            error('vdx only supports assigning lb, ub, or init for a variable via dot indexing');
                        end
                    else
                        error("Currently vdx allows only scalar assignment to lb, ub, or init via dot indexing")
                    end
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here.
                end
            end
            % TODO(@anton) dot indexng synatctic sugar needs some more thought and possibly a _lot_ more logic to handle
            %              different modalities of the RHS, including assigning to multiple indexes at once.
        end

        function n = parenListLength(obj,index_op,ctx)
           n = 1;
        end

        function varargout = braceReference(obj, index_op)
            if isscalar(index_op)
                values = cellfun(@(x) obj.solver.results.w(x), obj.indices, 'uni', false);
                out = squeeze(obj.symbolics.(index_op));
            else
                error('Brace indexing only accesses current values of the solution and no chained indexing')
            end
            varargout{1} = out;
        end

        function obj = braceAssign(obj,index_op,varargin)
            if isscalar(index_op)
            end
        end

        function n = braceListLength(obj,index_op,ctx)
           n = 1;
        end

        function obj = parenDelete(obj,index_op)
            error('Deletion of symbolics is not supported through the variable view')
        end

        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);

            cp.vector = [];
            cp.solver = [];
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            obj = NosnocVariable([],[]);
        end
    end
end
