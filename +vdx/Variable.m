classdef Variable < handle &...
        matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.indexing.RedefinesDot &...
        matlab.mixin.CustomDisplay &...
        matlab.mixin.Copyable
    properties
        % Indices of this :class:`vdx.Variable` in its :class:vdx.Vector.
        %
        %:type: cell
        indices

        % :class:`vdx.Vector` that this :class:`vdx.Variable` is a member of.
        vector

        % name of this variable
        name
    end

    properties (Hidden, SetAccess=private)
        reorder_to_end = false
    end

    properties (SetAccess=private)
        % Number of indices supported by this variable.
        %
        %:type: double
        depth = []
    end

    methods(Access=public)
        function obj = Variable(vector, name)
            obj.indices = cell(0,1);
            obj.vector = vector;
            obj.name = name;
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

        function output = to_string(obj, varargin)
        % Pretty prints this variable with the specified columns.
        %
        % Available columns are: 'sym', 'lb', 'ub', 'init', 'res', and 'mult', which are passed as string arguments to this method.
        % Default outputs all columns.

            printed_cols = [];
            % Calculate which cols to print.
            if isempty(varargin)
                printed_cols = [obj.vector.numerical_properties, obj.vector.numerical_outputs, "sym"];
            else
                printed_cols = [varargin{:}];
                if ismember("sym", printed_cols)
                    idx = find(printed_cols == "sym");
                    printed_cols(idx) = [];
                    printed_cols = [printed_cols "sym"];
                end
            end

            % Generate header
            header = '';
            header = [header sprintf('%-5s', 'i')];
            for name=printed_cols
                header = [header, sprintf('| %-12s', name)];
            end
            header = [header '\n'];

            % iterate over all requested values
            indices = sort([obj.indices{:}]);
            output = [];
            for ii=indices
                pline = sprintf('%-5d', ii);
                for name=printed_cols
                    if strcmp(name,"sym")
                        vec = obj.vector.sym;
                        pline = [pline '| ' char(formattedDisplayText(vec(ii)))];
                    else
                        vec = obj.vector.numerical_vectors.(name);
                        pline = [pline sprintf('| %-12.5g', vec(ii))];
                    end
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            output = [header output];
        end
        
        function dummy = print(obj, varargin)
        % Pretty prints this variable with the specified columns.
        %
        % Available columns are: 'sym', 'lb', 'ub', 'init', 'res', and 'mult', which are passed as string arguments to this method.
        % Default prints all columns.
            output = obj.to_string(varargin{:});
            fprintf(output);
            dummy = [];
        end

        function json = jsonencode(obj, varargin)
            obj.vector.apply_queued_assignments();
            var_struct = struct();

            var_struct.name = obj.name;
            var_struct.indices = permute(obj.indices, ndims(obj.indices):-1:1);
            var_struct.ind_shape = size(obj.indices);
            var_struct.depth = obj.depth;
            var_struct.reorder_to_end = obj.reorder_to_end;
            
            json = jsonencode(var_struct, varargin{:});
        end
    end

    methods (Static)
        function var = from_json(var_struct)
            var = vdx.Variable([],"");
            if ~iscell(var_struct.indices)
                var_struct.indices = {var_struct.indices};
            end
            var.indices = reshape(var_struct.indices, var_struct.ind_shape');
            var.name = var_struct.name;
            var.depth = var_struct.depth;
            var.reorder_to_end = var_struct.reorder_to_end;
        end
    end
    
    methods (Access=protected)
        function varargout = dotReference(obj,index_op)
            obj.vector.apply_queued_assignments();
            name = index_op(1).Name;
            if ~ismember(name,[obj.vector.numerical_properties, obj.vector.numerical_outputs])
                error(['numerical vector ' char(name) ' is does not exist for this vector'])
            end
            if ~isscalar(index_op)
                % TODO(@anton) better error here
                error('Compound indexing is not supported')
            end
            vec = obj.vector.(name);
            out = cellfun(@(x) vec(x), obj.indices, 'uni', false);
            out = permute(out, ndims(out):-1:1);
            varargout{1} = [out{:}];
        end

        function obj = dotAssign(obj,index_op,varargin)
            error('Assigning to unindexed variable is not supported')
        end

        function n = dotListLength(obj,index_op,indexContext)
            n=1;
        end
        
        function varargout = parenReference(obj, index_op)
            obj.vector.apply_queued_assignments();
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
                symbolics = cellfun(@(x) obj.vector.sym(x), obj.indices(adj_ind{:}), 'uni', false);
                out = squeeze(symbolics);
            else
                if index_op(2).Type == 'Dot'
                    if ~ismember(index_op(2).Name, [obj.vector.numerical_properties, obj.vector.numerical_outputs])
                        % TODO(@anton) better error
                        error('Attempt to access unsupported property for this vector type')
                    end
                    num_vec = obj.vector.numerical_vectors.(index_op(2).Name);
                    num = cellfun(@(x) num_vec(x), obj.indices(adj_ind{:}), 'uni', false);
                    out = squeeze(num);
                else
                    error('unsupported indexing');
                    % TODO(@anton) better error here.
                end
            end
            if isscalar(out)
                out = out{1};
            else
                out = permute(out, ndims(out):-1:1);
                out = [out{:}];
            end
            if numel(out) == 0
                out = [];
            end
            varargout{1} = out;
        end

        function obj = parenAssign(obj,index_op,varargin)
            if isempty(obj.depth)
                obj.depth = length(index_op(1).Indices);
                if obj.depth == 0
                    size_args = {0,1};
                else
                    size_args = num2cell([zeros(obj.depth, 1);1]);
                end
                obj.indices = cell(size_args{:});
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
                    indices = obj.vector.add_variable(index_op(1).Indices, arg{:});
                    adj_ind = index_adjustment(index_op.Indices);
                    obj.indices{adj_ind{:},1} = indices;
                elseif is_index_logical_array(index_op(1).Indices) % A boolean array representing where variables should be created
                    error('indexing via logical array not yet supported')
                else % Assume we want to assign multiple values
                    if false && iscell(arg{1}) && ischar(arg{1}{1}) % Fast variable creation sacrificing naming
                        inorderlst = all_combinations(index_op.Indices{:});
                        n_args = size(inorderlst, 1);
                        n_var = arg{1}{2};
                        arg{1}{2} = arg{1}{2}*n_args;
                        all_indices = obj.vector.add_variable_fast(n_args, arg{:});
                        
                        % create vars and assign.
                        for ii=1:n_args
                            curr = inorderlst(ii,:);
                            curr_cell = num2cell(curr);
                            adj_ind = index_adjustment(curr_cell);
                            indices = all_indices(1+(((ii-1)*n_var):(ii*n_var-1)));
                            obj.indices{adj_ind{:},1} = indices;
                        end
                    else 
                        inorderlst = all_combinations(index_op.Indices{:});

                        % create vars and assign.
                        for ii=1:size(inorderlst, 1)
                            curr = inorderlst(ii,:);
                            curr_cell = num2cell(curr);
                            adj_ind = index_adjustment(curr_cell);
                            indices = obj.vector.add_variable(curr, arg{:});
                            adj_ind = index_adjustment(curr_cell);
                            obj.indices{adj_ind{:},1} = indices;
                        end
                    end
                end
            else
                if index_op(2).Type == 'Dot'
                    obj.vector.apply_queued_assignments();
                    if is_index_scalar(index_op(1).Indices)
                        if ~ismember(index_op(2).Name, [obj.vector.numerical_properties])
                            % TODO(@anton) better error
                            error('Attempt to assign to unsupported property for this vector type')
                        end
                        adj_ind = index_adjustment(index_op(1).Indices);
                        obj.vector.numerical_vectors.(index_op(2).Name)(obj.indices{adj_ind{:}}) = varargin{1};
                    else
                        % TODO (@anton) preempt wrong size with a good error
                        % TODO (@anton) allow for structured right hand side
                        inorderlst = all_combinations(index_op(1).Indices{:});
                        n_args = size(inorderlst, 1);
                        for ii=1:n_args
                            curr = inorderlst(ii,:);
                            curr_cell = num2cell(curr);
                            adj_ind = index_adjustment(curr_cell);
                            obj.vector.numerical_vectors.(index_op(2).Name)(obj.indices{adj_ind{:}}) = varargin{1};
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

        function propgrp = getPropertyGroups(obj)
            gTitle1 = 'Numeric properties';
            gTitle2 = 'Numeric Outputs';
            propList1 = struct;
            for name=obj.vector.numerical_properties
                propList1.(name) = [];% TODO(@anton) better custom display here
            end
            propList2 = struct;
            for name=obj.vector.numerical_outputs
                propList2.(name) = [];% TODO(@anton) better custom display here
            end
            propgrp(1) = matlab.mixin.util.PropertyGroup(propList1,gTitle1);
            propgrp(2) = matlab.mixin.util.PropertyGroup(propList2,gTitle2);
        end

        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            scalarHeader = [className];
            header = sprintf('%s\n',scalarHeader);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end

        function displayNonScalarObject(obj)
            dimStr = matlab.mixin.CustomDisplay.convertDimensionsToString(obj);
            cName = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            headerStr = [dimStr,' ',cName];
            header = sprintf('%s\n',headerStr);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end
    end

    methods
        function is_terminal = set_is_terminal(obj, is_terminal)
        % This only does anything for scalar variables
            obj.reorder_to_end = is_terminal;
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
