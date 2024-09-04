classdef Vector < handle &...
        matlab.mixin.indexing.RedefinesDot &...
        matlab.mixin.CustomDisplay &...
        dynamicprops &...
        matlab.mixin.Copyable
% A class which provides a wrapper around CasADi symbolics and tracks the indicies of :class:`vdx.Variable` within it.
%
% :param vdx.Problem problem: Problem which this vector is a member of.
% :param string casadi_type: either 'SX' (default) or 'MX' which determines the kind of CasADi symbolic stored.
    properties (Access=public)
        % CasADi symbolic vector that is wrapped by this object.
        %
        %:type: casadi.SX|casadi.MX
        sym
        
        % Casadi type
        casadi_type

        % pointer to parent problem
        problem
    end
    
    properties (Access=protected)
        % Internal struct of index tracking variables
        variables struct

        % Assignment queue
        pending_assignments

        % current length with pending assignements
        len
    end

    properties (Access={?vdx.Variable,?vdx.Vector})
         % internal struct of numerical data associated with this
        numerical_vectors struct
    end
    
    properties (Access=protected, Abstract)
        % default values for numerical vectors
        numerical_defaults
    end

    properties (Constant, Abstract, Hidden)
        numerical_properties
        numerical_outputs
        allow_nonscalar_symbolics
        allow_nonsymbolic_assignment
    end

    methods (Access=public)
        function obj = Vector(problem, varargin)
            p = inputParser;
            addRequired(p, 'problem');
            for name=obj.numerical_properties
                if isfield(obj.numerical_defaults, name)
                    default = obj.numerical_defaults.(name);
                else
                    default = 0;
                end
                addParameter(p, ['default_' char(name)], default);
            end
            addParameter(p, 'casadi_type', 'SX');
            parse(p, problem, varargin{:});

            % Populate core functionality
            obj.sym = casadi.(p.Results.casadi_type);
            obj.problem = problem;
            obj.variables = struct;
            obj.numerical_vectors = struct;
            obj.casadi_type = p.Results.casadi_type;
            obj.len = 0;

            % Populate defaults
            for name=obj.numerical_properties
                obj.numerical_defaults.(name) = p.Results.(['default_' char(name)]);
            end

            % Populate vectors
            for name=obj.numerical_properties
                obj.numerical_vectors.(name) = [];
            end
            for name=obj.numerical_outputs
                obj.numerical_vectors.(name) = [];
            end
        end

        function print(obj, varargin)
        % Pretty prints this vector with the specified columns.
        %
        % Available columns are the strings in the union of :attr:`numerical_properties` and :attr:`numerical_outputs`, which are passed as string arguments to this method.
        % Default prints all columns.
            output = obj.to_string(varargin{:});
            fprintf(output);
        end

        function output = to_string(obj, varargin)
            obj.apply_queued_assignments();
            printed_cols = [];
            % Calculate which cols to print.
            if isempty(varargin)
                printed_cols = [obj.numerical_properties, obj.numerical_outputs, "sym"];
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
            n = size(obj.sym, 1);
            output = '';
            for ii=1:n
                pline = sprintf('%-5d', ii);
                for name=printed_cols
                    if strcmp(name,"sym")
                        pline = [pline '| ' char(formattedDisplayText(obj.sym(ii)))];
                    else
                        vec = obj.numerical_vectors.(name);
                        pline = [pline sprintf('| %-12.5g', vec(ii))];
                    end
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            output = [header output];
        end

        function order_indices = sort_by_index(obj)
        % Sorts this vector so that the vectors occur in column major order with lower dimensional :class:`vdx.Variable`
        % always occuring before higher dimensional :class:`vdx.Variable`. This can be useful to recover any sparsity structure in the 
        % problem that comes from the structure of constaraints and variables.
            obj.apply_queued_assignments();
            vars = fieldnames(obj.variables);
            % get depth
            lengths = 1;
            for ii=1:numel(vars)
                if isa(obj.variables.(vars{ii}), 'vdx.VariableGroup')
                    continue
                end
                s = size(obj.variables.(vars{ii}));
                dims = ndims(obj.variables.(vars{ii}));
                ls = length(s);
                ll = length(lengths);
                if ls > ll
                    lengths = max([lengths, zeros(1,ls-ll)], s);
                elseif ls < ll
                    lengths = max(lengths, [s, zeros(1,ll-ls)]);
                else
                    lengths = max(lengths, s);
                end
            end

            vars_by_depth = cell(length(lengths)+1, 1);
            for ii=1:length(vars_by_depth)
                vars_by_depth{ii} = {};
            end
            for ii=1:numel(vars)
                if isa(obj.variables.(vars{ii}), 'vdx.VariableGroup')
                    continue
                end
                dims = obj.variables.(vars{ii}).depth+1;
                vars_by_depth{dims} = vertcat(vars_by_depth{dims}, vars(ii));
            end
            
            indices = {};
            for len=lengths
                indices = horzcat(indices, {1:len});
            end
            % build modlist
            modlist = cumprod([1,flip(lengths)]);
            modlist = flip(modlist(1:end-1));
            
            % We subtract 1 to get the 0 indexing correct :)
            inorderlst = all_combinations(indices{:})-1;

            order_indices = zeros(size(obj.sym));
            n_new = 0;
            
            % First re-normalize 0 dimensional vars (i.e indicies 1x1)
            d_vars = vars_by_depth{1};
            terminal_vars = {};
            for jj=1:numel(d_vars)
                var = obj.variables.(d_vars{jj});
                if var.reorder_to_end
                    terminal_vars = [terminal_vars, {var}];
                    continue;
                end
                
                ind = var.indices{1};
                n = numel(ind);
                indices = (n_new+1):(n_new+n);
                n_new = n_new + n;
                order_indices(indices) = ind;
                
                var.indices{1} = indices;
            end
            
            % now handle rest of vars
            dims = 1:length(lengths);
            for ii=1:size(inorderlst,1)
                % This is an ugly hack. we create a cell array in order to create a comma separated list for indexing
                curr = inorderlst(ii,:);
                mods = mod(ii-1,modlist);
                mask = mods == 0;
                dims_to_process = dims(mask);
                for dim=dims_to_process
                    d_vars = vars_by_depth{dim+1};
                    curr_for_dim = num2cell(curr(1:dim));
                    for jj=1:numel(d_vars)
                        var = obj.variables.(d_vars{jj});
                        curr_for_dim_adj = num2cell(curr(1:dim)+1);

                        try
                            ind = var.indices{curr_for_dim_adj{:}};
                            n = numel(ind);
                            indices = (n_new+1):(n_new+n);
                            n_new = n_new + n;
                            order_indices(indices) = ind;
                            
                            var.indices{curr_for_dim_adj{:}} = indices;
                        end
                    end
                end
            end

            % put terminal vars at the end
            for jj=1:numel(terminal_vars)
                var = terminal_vars{jj};
                
                ind = var.indices{1};
                n = numel(ind);
                indices = (n_new+1):(n_new+n);
                n_new = n_new + n;
                order_indices(indices) = ind;
                
                var.indices{1} = indices;
            end
            
            obj.sym = obj.sym(order_indices);
            for name=obj.numerical_properties
                obj.numerical_vectors.(name) = obj.numerical_vectors.(name)(order_indices);
            end
        end

        function add_variable_group(obj, name, vars, varargin)
        % Adds a :class:`vdx.VariableGroup` to this vector
            if isfield(obj.variables,name)
                error('Variable or VariableGroup with this name already exists')
            else
                obj.variables.(name) = vdx.VariableGroup(vars, varargin{:});
            end
        end

        function has = has_var(obj, name)
            has = isfield(obj.variables, char(name));
        end

        function apply_queued_assignments(obj)
            if isempty(obj.pending_assignments)
                return
            end
            obj.sym = vertcat(obj.sym, obj.pending_assignments.syms);

            for name=obj.numerical_properties
                obj.numerical_vectors.(name) = vertcat(obj.numerical_vectors.(name), obj.pending_assignments.(name));
            end

            for name=obj.numerical_outputs
                obj.numerical_vectors.(name) = vertcat(obj.numerical_vectors.(name), obj.pending_assignments.(name));
            end

            obj.pending_assignments = [];
        end

        function vars = get_vars(obj)
            vars = struct2cell(obj.variables);
        end

        function vars = get_var_names(obj)
            vars = fieldnames(obj.variables);
        end
    end
    
    methods (Access=protected)
        function varargout = dotReference(obj,index_op)
            obj.apply_queued_assignments();
            name = index_op(1).Name;
            % first try numerical vectors
            if ismember(name,[obj.numerical_properties, obj.numerical_outputs])
                varargout{1} = obj.numerical_vectors.(name);
                return
            end
            if ~isfield(obj.variables, name)
                err.message = sprintf(['Variable ' char(name) ' does not exist on this vector']);
                err.identifier = 'vdx:indexing:assign_to_scalar';
                stack = dbstack('-completenames');
                stack(1).name = 'Vector.reference';
                err.stack = stack(1:end);
                error(err);
                error([])
            end
            if isscalar(index_op)
                warning('on', 'verbose')
                warning('on', 'backtrace')
                warning('vdx:indexing:dot_reference_returns_vdx_var',...
                    sprintf(['You have accessed vdx.Variable ' char(name) ' directly.\n'...
                              'In many cases this is a mistake unless you are using advanced indexing features of vdx.']));
                warning('off', 'backtrace')
                warning('off', 'verbose')
            end
            varargout{1} = obj.variables.(index_op);
        end

        function obj = dotAssign(obj,index_op,varargin)
            if string(version('-release')) < "2023a"
                persistent matlab_bug_workaround_to_prevent_garbage_collection;
                matlab_bug_workaround_to_prevent_garbage_collection = varargin{1};
            end
            name = index_op(1).Name;
            % first try numerical vectors
            if ismember(name,[obj.numerical_properties, obj.numerical_outputs])
                obj.apply_queued_assignments();
                obj.numerical_vectors.(index_op) = varargin{1};
                return
            end
            if isscalar(index_op)
                if ~isfield(obj.variables,name) % Workaround for scalar variables because matlab throws a fit if you try x() = 1;
                    var = vdx.Variable(obj, name);
                    obj.variables.(name) = var;
                    P = obj.addprop(name);
                    obj.(name) = var;
                    var(vdx.constants.scalar{:}) = varargin{1};
                    return
                elseif obj.variables.(indexop(1)).depth == 0; % TODO maybe this should also error.
                    var = obj.variables.(indexop(1));
                    var(vdx.constants.scalar{:}) = varargin{1};
                else
                    err.message = sprintf(['Assigning directly to variable ' char(name) ' is not allowed. Include an index.\n'...
                                            'If you want to track a scalar variable (with no subscripts you can index via: ' char(name) '()']);
                    err.identifier = 'vdx:indexing:assign_to_scalar';
                    stack = dbstack('-completenames');
                    stack(1).name = 'Vector.assignment';
                    err.stack = stack;
                    error(err);
                end
            end
            if ~isfield(obj.variables,name)
                var = vdx.Variable(obj, name);
                obj.variables.(name) = var;
                P = obj.addprop(name);
                obj.(name) = var;
            end
            obj.variables.(index_op) = varargin{1};
        end
        
        function n = dotListLength(obj,index_op,indexContext)
            n=1;
        end

        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);

            % deepcopy variables
            var_names = fieldnames(obj.variables);
            for ii=1:length(var_names)
                cp.variables.(var_names{ii}) = copy(obj.variables.(var_names{ii}));
                cp.variables.(var_names{ii}).vector = cp;
            end
        end

        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                gTitle1 = 'Numeric properties';
                gTitle2 = 'Numeric Outputs';
                propList1 = struct;
                for name=obj.numerical_properties
                    propList1.(name) = obj.numerical_vectors.(name);
                end
                propList2 = struct;
                for name=obj.numerical_outputs
                    propList2.(name) = obj.numerical_vectors.(name);
                end
                propgrp(1) = matlab.mixin.util.PropertyGroup({'sym'});
                propgrp(2) = matlab.mixin.util.PropertyGroup(propList1,gTitle1);
                propgrp(3) = matlab.mixin.util.PropertyGroup(propList2,gTitle2);
                propgrp(4) = matlab.mixin.util.PropertyGroup(obj.variables, 'Variables');
            end
        end

        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            scalarHeader = [className];
            header = sprintf('%s\n',scalarHeader);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end
    end

    methods (Access={?vdx.Variable, ?vdx.VariableGroup, ?vdx.ConstraintVector})
        function indices = add_variable(obj, indices, symbolic, varargin)
        % Adds a :class:`vdx.Variable` to the internal symbolic and numeric vectors.
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            for name=obj.numerical_properties
                default = obj.numerical_defaults.(name);
                addOptional(p, name, default);
            end
            parse(p, obj, symbolic, varargin{:});

            symbolic = p.Results.symbolic;
            numeric_vals = struct;
            for name=obj.numerical_properties
                astruct.(name) = p.Results.(name);
            end

            % Check that symbolic is valid
            if iscell(symbolic) && ~isa(symbolic{1}, 'casadi.Function') && ~obj.allow_nonsymbolic_assignment
                % TODO better error
                error(['This vector of class ' class(obj) ' does not allow for {name, size} form of assignment.'])
            end
            % TODO more checks here
            
            % Handle non-symbolic input as (name, size) pair
            symbolic = obj.eval_symbolic(symbolic, indices);

            % Get size and populate possibly empty values
            n = size(symbolic, 1);
            for name=obj.numerical_properties
                if length(p.Results.(name)) == 1
                    astruct.(name) = astruct.(name) * ones(n,1);
                end
            end
            
            lens = [size(symbolic,1)];
            for name=obj.numerical_properties
                lens = [lens, size(astruct.(name), 1)];
            end
            if ~all(lens == lens(1))
                % TODO(@anton) better error message
                error("mismatched dims")
            end

            for name=obj.numerical_outputs
                astruct.(name) = zeros(n,1);
            end

            n_sym = obj.len;

            indices = (n_sym+1):(n_sym+n);
            obj.len = obj.len + n;
            astruct.syms = symbolic;
            if isempty(obj.pending_assignments)
                obj.pending_assignments = astruct;
            else
                obj.pending_assignments(end+1) = astruct;
            end
        end

        function indices = add_variable_fast(obj, n_args, symbolic, varargin)
        % Adds a :class:`vdx.Variable` to the internal symbolic and numeric vectors.
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            for name=obj.numerical_properties
                default = obj.numerical_defaults.(name);
                addOptional(p, name, default);
            end
            parse(p, obj, symbolic, varargin{:});

            symbolic = p.Results.symbolic;
            numeric_vals = struct;
            for name=obj.numerical_properties
                numeric_vals.(name) = p.Results.(name);
            end

            % Check that symbolic is valid
            if iscell(symbolic) && ~isa(symbolic{1}, 'casadi.Function') && ~obj.allow_nonsymbolic_assignment
                % TODO better error
                error(['This vector of class ' class(obj) ' does not allow for {name, size} form of assignment.'])
            end
            % TODO more checks here
            
            % Handle non-symbolic input as (name, size) pair
            symbolic = obj.eval_symbolic(symbolic);

            % Get size and populate possibly empty values
            n = size(symbolic, 1);
            for name=obj.numerical_properties
                if length(p.Results.(name)) == 1
                    numeric_vals.(name) = numeric_vals.(name) * ones(n,1);
                else
                    numeric_vals.(name) = repmat(numeric_vals.(name), n_args, 1);
                end
            end
            
            lens = [size(symbolic,1)];
            for name=obj.numerical_properties
                lens = [lens, size(numeric_vals.(name), 1)];
            end
            if ~all(lens == lens(1))
                % TODO(@anton) better error message
                error("mismatched dims")
            end

            n_sym = size(obj.sym, 1);

            obj.sym = vertcat(obj.sym, symbolic);
            for name=obj.numerical_properties
                obj.numerical_vectors.(name) = [obj.numerical_vectors.(name); numeric_vals.(name)];
            end

            % initialize results and multipliers to zero
            % TODO(@anton) is there a better descision than this?
            for name=obj.numerical_outputs
                numerics.(name) = zeros(n,1);
                obj.numerical_vectors.(name) = [obj.numerical_vectors.(name); zeros(n,1)];
            end
            indices = (n_sym+1):(n_sym+n);
        end
    end

    methods (Access=private)
        function sym = eval_symbolic(obj, sym, varargin)
        % process sym into a casadi symbolic possibly renaming using indices in the process
            p = inputParser;
            addRequired(p, 'sym');
            addOptional(p, 'index', []);
            parse(p, sym, varargin{:});

            % Pile of ifs to handle different cases
            if isempty(p.Results.index)
                if any(size(sym) == 0) || isa(sym, ['casadi.' obj.casadi_type])
                    % do nothing, directly pass through
                elseif iscell(sym) && length(sym) == 2 &&...
                        ischar(sym{1}) && length(sym{2}) == 1 && isnumeric(sym{2}) && round(sym{2}) == sym{2}
                    sym = define_casadi_symbolic(obj.casadi_type, sym{1}, sym{2});
                else
                    error("Incorrect type")
                end
            else
                if any(size(sym) == 0)
                    % Do nothing
                elseif isa(sym, ['casadi.' obj.casadi_type])
                    if sym.is_symbolic
                        % if we can make the names of symbolics nice. However this may be bad if one wants to add single variable constraints in g instead of w.
                        % This is bugged do nothing instead
                        % name = split(sym(1).name, '_');
                        % name = [name{1:end-1} index_string(p.Results.index)];
                        % sym = define_casadi_symbolic(obj.casadi_type, name, size(sym, 1));
                    else
                        % pass through and do nothing
                    end
                elseif iscell(sym) && length(sym) >= 2
                    if ischar(sym{1}) && length(sym{2}) == 1 && round(sym{2}) == sym{2}
                        name = [sym{1} index_string(p.Results.index)];
                        sym = define_casadi_symbolic(obj.casadi_type, name, sym{2});
                    elseif isa(sym{1}, 'casadi.Function')
                        varargs = sym(3:end);
                        arg_group = vdx.VariableGroup(sym{2}, varargs{:});
                        fun = sym{1};
                        inds = num2cell(p.Results.index);
                        sym_args = arg_group{inds{:}};
                        sym = fun(sym_args{:});
                    else
                        error("Incorrect type")
                    end 
                else
                    error("Incorrect type")
                end
            end
        end
    end
end

function res = is_index_scalar(index)
    res = all(cellfun(@(x) isscalar(x) & ~ischar(x), index));
end

function sym = define_casadi_symbolic(type, name, dims, sparsity)
    import casadi.*
    if nargin < 3
        dims = 1;
    end

    if ~exist('sparsity', 'var')
        if size(dims, 2) == 1
            sparsity = Sparsity.dense([dims,1]);
        else
            sparsity = Sparsity.dense(dims);
        end
    end

    if strcmp(type, 'casadi.SX') || strcmp(type, 'SX')
        sym = SX.sym(name, sparsity);
    elseif strcmp(type, 'casadi.MX')  || strcmp(type, 'MX')
        sym = MX.sym(name, sparsity);
    else
        error('Type must be MX or SX.')
    end
end

function valid = valid_index(var,idx)
    ni = length(idx);
    nv = ndims(idx);
end

function istring = index_string(indices)
    istring = '';
    if iscell(indices)
        indices = [indices{:}];
    end
    for i=indices
        istring = [istring '_' num2str(i)];
    end
end
