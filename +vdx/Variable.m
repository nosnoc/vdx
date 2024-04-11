classdef Variable < handle &...
        matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.Copyable
    properties
        % Indices of this :class:`vdx.Variable` in its :class:vdx.Vector.
        %
        %:type: cell
        indices

        % :class:`vdx.Vector` that this :class:`vdx.Variable` is a member of.
        vector
    end

    properties (Dependent)
        % Horizonally concatenated vectors of all lower bounds of this :class:`vdx.Variable`.
        %
        %:type: double
        lb

        % Horizonally concatenated vectors of all upper bounds of this :class:`vdx.Variable`.
        %
        %:type: double
        ub

        % Horizonally concatenated vectors of all initial values of this :class:`vdx.Variable`.
        %
        %:type: double
        init

        % Horizonally concatenated vectors of all results of this :class:`vdx.Variable`.
        %
        %:type: double
        res

        % Horizonally concatenated vectors of all Lagrange multipliers of this :class:`vdx.Variable`.
        %
        %:type: double
        mult
    end

    properties (SetAccess=private)
        % Number of indices supported by this variable.
        %
        %:type: double
        depth = []
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
            if obj.depth ~= length(index_op(1).Indices)
                err.message = sprintf(['You are subscripting a variable using ' num2str(length(index_op(1).Indices)) ' subscripts but this variable expects ' num2str(obj.depth) ' subscripts.']);
                err.identifier = 'vdx:indexing:incorrect_num_of_subscripts';
                stack = dbstack('-completenames');
                stack(1).name = 'Variable.reference';
                err.stack = stack;
                error(err);
            end
            if ~all_indices_integral(index_op(1).Indices)
                err.message = sprintf(['Variable subscripts must be integral.']);
                err.identifier = 'vdx:indexing:non_integral_subscripts';
                stack = dbstack('-completenames');
                stack(1).name = 'Variable.reference';
                err.stack = stack;
                error(err);
            end
            adj_ind = vdx.indexing.identity(index_op(1).Indices, obj);
            unacceptable_indices = obj.unacceptable_indices(adj_ind);
            if ~isempty(unacceptable_indices)
                err.message = sprintf(['Subscript in position(s) [' num2str(unacceptable_indices) '] exceed(s) maximum subscripts']);
                err.identifier = 'vdx:indexing:non_integral_subscripts';
                stack = dbstack('-completenames');
                stack(1).name = 'Variable.reference';
                err.stack = stack;
                error(err);
            end
            if isscalar(index_op)
                symbolics = cellfun(@(x) obj.vector.sym(x), obj.indices, 'uni', false);
                out = squeeze(symbolics(adj_ind{:}));
            else
                if index_op(2).Type == 'Dot'
                    switch(index_op(2).Name)
                      case "lb"
                        lb = cellfun(@(x) obj.vector.lb(x), obj.indices, 'uni', false);
                        out = squeeze(lb(adj_ind{:}));
                      case "ub"
                        ub = cellfun(@(x) obj.vector.ub(x), obj.indices, 'uni', false);
                        out = squeeze(ub(adj_ind{:}));
                      case "init"
                        init = cellfun(@(x) obj.vector.init(x), obj.indices, 'uni', false);
                        out = squeeze(init(adj_ind{:}));
                      case "res"
                        res = cellfun(@(x) obj.vector.res(x), obj.indices, 'uni', false);
                        out = squeeze(res(adj_ind{:}));
                      case "mult"
                        mult = cellfun(@(x) obj.vector.mult(x), obj.indices, 'uni', false);
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
            if isempty(obj.depth)
                obj.depth = length(index_op(1).Indices);
            end
            if obj.depth ~= length(index_op(1).Indices)
                err.message = sprintf(['You are assigning to variable using ' num2str(length(index_op(1).Indices)) ' subscripts but this variable expects ' num2str(obj.depth) ' subscripts.']);
                err.identifier = 'vdx:indexing:incorrect_num_of_subscripts';
                stack = dbstack('-completenames');
                stack(1).name = 'Variable.assignment';
                err.stack = stack;
                error(err);
            end
            if isscalar(index_op)
                % get the cell array of args (x,x0,lbx,ubx)
                arg = varargin{1};
                % allow for multi index variable creation
                if is_index_scalar(index_op(1).Indices) % A single scalar variable     
                    indices = obj.vector.add_variable([index_op(1).Indices{:}], arg{:});
                    adj_ind = index_adjustment(index_op.Indices);
                    obj.indices{adj_ind{:},1} = indices;
                elseif is_index_logical_array(index_op(1).Indices) % A boolean array representing where variables should be created
                    error('indexing via logical array not yet supported')
                else % Assume we want to assign multiple values
                    inorderlst = all_combinations(index_op.Indices{:});

                    % create vars and assign.
                    for ii=1:size(inorderlst)
                        curr = inorderlst(ii,:);
                        curr_cell = num2cell(curr);
                        adj_ind = index_adjustment(curr_cell);
                        indices = obj.vector.add_variable(curr, arg{:});
                        adj_ind = index_adjustment(curr_cell);
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

    methods (Access=private)
        function unacceptable = unacceptable_indices(obj, adj_ind)
            sz = size(obj.indices);
            sz = sz(1:obj.depth);
            if is_index_scalar(adj_ind)
                unacceptable = find([adj_ind{:}] > sz);
            else
                more_adj_ind = adj_ind;
                for ii=1:length(adj_ind)
                    if strcmp(adj_ind{ii}, ':')
                        more_adj_ind{ii} = 1;
                    end
                end
                inorderlst = all_combinations(more_adj_ind{:});

                unacceptable = find(sum(inorderlst > repmat(sz, size(inorderlst,1), 1)), 1);
            end
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

function res = all_indices_integral(indices)
    res = all(cellfun(@(x) (isnumeric(x) & all(round(x) == x)) | strcmp(x,':'), indices));
end
