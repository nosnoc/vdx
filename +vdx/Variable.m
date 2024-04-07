classdef Variable < handle &...
        matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.Copyable
    properties
        % Indices of this :class: vdx.Variable in its :class: vdx.Vector.
        %
        %:type: cell
        indices

        % :class: vdx.Vector that this :class: vdx.Variable is a member of.
        vector
    end

    properties (Dependent)
        % Horizonally concatenated vectors of all lower bounds of this :class: vdx.Variable.
        %
        %:type: double
        lb

        % Horizonally concatenated vectors of all upper bounds of this :class: vdx.Variable.
        %
        %:type: double
        ub

        % Horizonally concatenated vectors of all initial values of this :class: vdx.Variable.
        %
        %:type: double
        init

        % Horizonally concatenated vectors of all results of this :class: vdx.Variable.
        %
        %:type: double
        res

        % Horizonally concatenated vectors of all Lagrange multipliers of this :class: vdx.Variable.
        %
        %:type: double
        mult

        % Number of indices supported by this variable.
        %
        %:type: double
        depth
    end
    
    methods
        function out = get.lb(obj)
            out = cellfun(@(x) obj.vector.lb(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            out = [out{:}];
        end
        function out = get.ub(obj)
            out = cellfun(@(x) obj.vector.ub(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            out = [out{:}];
        end
        function out = get.init(obj)
            out = cellfun(@(x) obj.vector.init(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            out = [out{:}];
        end
        function out = get.res(obj)
            out = cellfun(@(x) obj.vector.res(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            out = [out{:}];
        end
        function out = get.mult(obj)
            out = cellfun(@(x) obj.vector.mult(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            out = [out{:}];
        end

        function out = get.depth(obj)
            sz = size(obj.indices);
            out = sum(sz > 1);
        end
    end

    methods(Access=public)
        function obj = Variable(vector)
            obj.indices = cell(0,1);
            obj.vector = vector;
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
        
        function output = print(obj, varargin)
        % Pretty prints this variable with the specified columns.
        %
        % Available columns are: 'sym', 'lb', 'ub', 'init', 'res', and 'mult', which are passed as string arguments to this method.
        % Default prints all columns.
            sym = false;
            lb = false;
            ub = false;
            init = false;
            res = false;
            mult = false;
            if isempty(varargin)
                sym = true;
                lb = true;
                ub = true;
                init = true;
                res = true;
                mult = true;
            else
                if any(ismember(lower(varargin), 'sym'))
                    sym = true;
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
            if sym
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
                if sym
                    pline = [pline char(formattedDisplayText(obj.vector.sym(ii)))];
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
                adj_ind = vdx.indexing.identity(index_op.Indices, obj);
                symbolics = cellfun(@(x) obj.vector.sym(x), obj.indices(adj_ind{:}), 'uni', false);
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
                % allow for multi index variable creation
                if is_index_scalar(index_op(1).Indices) % A single scalar variable     
                    symbolic = arg{1}; % TODO this living here implies we should move other handeling out of 'add_variable'
                    if iscell(symbolic) && isa(symbolic{1}, 'casadi.Function')
                        varargs = symbolic(3:end);
                        arg_group = vdx.VariableGroup(symbolic{2}, varargs{:});
                        fun = symbolic{1};
                        fargs = arg_group{index_op(1).Indices{:}};
                        symbolic = fun(fargs{:});
                        arg{1} = symbolic;
                    end
                    indices = obj.vector.add_variable(arg{:});
                    adj_ind = index_adjustment(index_op.Indices);
                    obj.indices{adj_ind{:},1} = indices;
                elseif is_index_logical_array(index_op(1).Indices) % A boolean array representing where variables should be created
                    error('indexing via logical array not yet supported')
                else % Assume we want to assign multiple values
                    inorderlst = all_combinations(index_op.Indices{:});
                    % Check first argument and adjust naming
                    x = arg{1};
                    is_fun = false;
                    if iscell(x)
                        if isa(x{1}, 'casadi.Function')
                            is_fun = true;
                            varargs = x(3:end);
                            arg_group = vdx.VariableGroup(x{2}, varargs{:});
                            fun = x{1};
                        else
                            name = x{1}; n = x{2};
                        end
                    else % assume it is an SX or MX
                        name = x.name; n = size(1, x);
                    end

                    % create vars and assign.
                    for ii=1:size(inorderlst)
                        curr = inorderlst(ii,:);
                        curr_cell = num2cell(curr);
                        adj_ind = index_adjustment(curr_cell);
                        if is_fun
                            fargs = arg_group{curr_cell{:}};
                            symbolic = fun(fargs{:});
                            arg{1} = symbolic;
                        else
                            arg{1} = {[name index_string(curr)], n};
                        end
                        % add variable
                        indices = obj.vector.add_variable(arg{:});
                        obj.indices{adj_ind{:},1} = indices;
                    end
                end
            else
                if index_op(2).Type == 'Dot'
                    if is_index_scalar(index_op(1).Indices)
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
                          case "mult"
                            adj_ind = index_adjustment(index_op(1).Indices);
                            obj.vector.mult(obj.indices{adj_ind{:}}) = varargin{1};
                          otherwise
                            error('vdx only supports assigning lb, ub, mult, or init for a variable via dot indexing');
                        end
                    else
                        % TODO (@anton) preempt wrong size with a good error
                        % TODO (@anton) allow for structured right hand side
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
                          case "mult"
                            adj_ind = index_adjustment(index_op(1).Indices);
                            obj.vector.mult(obj.indices{adj_ind{:}}) = varargin{1};
                          otherwise
                            error('vdx only supports assigning lb, ub, mult, or init for a variable via dot indexing');
                        end
                    end
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here, (maybe do some heuristic on what user was trying to do).
                end
            end
            % TODO(@anton) dot indexng synatctic sugar needs some more thought and possibly a _lot_ more logic to handle
            %              different modalities of the RHS, including assigning to multiple indexes at once.
        end

        function n = parenListLength(obj,index_op,ctx)
           n = 1;
        end

        function obj = parenDelete(obj,index_op)
            error('Deletion of symbolics is not supported through the variable view')
        end

        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);

            cp.vector = [];
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            obj = [];
        end
    end
end

function res = is_index_scalar(index)
    res = all(cellfun(@(x) isscalar(x) & ~ischar(x), index));
end

function res = is_index_logical_array(index)
    if length(index) == 1 && islogical(index{1})
        res = true;
    else
        res = false;
    end
end

function istring = index_string(indices)
    istring = '';
    for i=indices
        istring = [istring '_' num2str(i)];
    end
end
