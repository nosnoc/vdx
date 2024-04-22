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
    end

    properties (SetAccess=private)
        % Number of indices supported by this variable.
        %
        %:type: double
        depth = []
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

            printed_cols = [];
            % Calculate which cols to print.
            if isempty(varargin)
                printed_cols = [obj.vector.numerical_properties, obj.vector.numerical_outputs, "sym"];
            else
                printed_cols = [varargin{:}];
            end

            % Generate header
            header = 'i\t\t';
            for name=printed_cols
                header = [header, char(name), '\t\t'];
            end
            header = [header '\n'];
            fprintf(header);

            % iterate over all requested values
            indices = sort([obj.indices{:}]);
            output = [];
            for ii=indices
                pline = [num2str(ii) '\t\t'];
                for name=printed_cols
                    if strcmp(name,"sym")
                        vec = obj.vector.sym;
                        pline = [pline char(formattedDisplayText(vec(ii)))];
                    else
                        vec = obj.vector.numerical_vectors.(name);
                        pline = [pline sprintf('%-8.5g\t', vec(ii))];
                    end
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            fprintf(output);
        end
    end

    methods (Access=protected)
        function varargout = dotReference(obj,index_op)
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
                    num = cellfun(@(x) num_vec(x), obj.indices, 'uni', false);
                    out = squeeze(num(adj_ind{:}));
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
                    indices = obj.vector.add_variable([], arg{:});
                    adj_ind = index_adjustment(index_op.Indices);
                    obj.indices{adj_ind{:},1} = indices;
                elseif is_index_logical_array(index_op(1).Indices) % A boolean array representing where variables should be created
                    error('indexing via logical array not yet supported')
                else % Assume we want to assign multiple values
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
            else
                if index_op(2).Type == 'Dot'
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
                        adj_ind = index_adjustment(index_op(1).Indices);
                        obj.vector.numerical_vectors.(index_op(2).Name)(obj.indices{adj_ind{:}}) = varargin{1};
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
