classdef Vector < handle &...
        matlab.mixin.indexing.RedefinesDot &...
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
    end
    
    properties (Access=private)
        % pointer to parent problem
        problem

        % Internal struct of index tracking variables
        variables struct

        % Casadi type
        casadi_type
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
            printed_cols = [];
            % Calculate which cols to print.
            if isempty(varargin)
                printed_cols = [obj.numerical_properties, obj.numerical_outputs, "sym"];
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
            n = size(obj.sym, 1);
            output = [];
            for ii=1:n
                pline = [num2str(ii) '\t\t'];
                for name=printed_cols
                    if strcmp(name,"sym")
                        pline = [pline char(formattedDisplayText(obj.sym(ii)))];
                    else
                        vec = obj.numerical_vectors.(name);
                        pline = [pline sprintf('%-8.5g\t', vec(ii))];
                    end
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            fprintf(output);
        end

        function sort_by_index(obj)
        % Sorts this vector so that the vectors occur in column major order with lower dimensional :class:`vdx.Variable`
        % always occuring before higher dimensional :class:`vdx.Variable`. This can be useful to recover any sparsity structure in the 
        % problem that comes from the structure of constaraints and variables.
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

            % new vectors.
            new_sym = [];
            new_numerics = struct;
            for name=[obj.numerical_properties, obj.numerical_outputs]
                new_numerics.(name) = [];
            end

            % First re-normalize 0 dimensional vars (i.e indicies 1x1)
            d_vars = vars_by_depth{1};
            for jj=1:numel(d_vars)
                var = obj.variables.(d_vars{jj});
                n_new_sym = length(new_sym);
                v_sym = var();
                new_sym = vertcat(new_sym, v_sym);
                for name=[obj.numerical_properties, obj.numerical_outputs]
                    new_numerics.(name) = vertcat(new_numerics.(name), var().(name));
                end

                n = length(v_sym);
                indices = (n_new_sym+1):(n_new_sym+n);
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
                        n_new_sym = length(new_sym);
                        v_sym = var(curr_for_dim{:});
                        new_sym = vertcat(new_sym, v_sym);
                        for name=[obj.numerical_properties, obj.numerical_outputs]
                            new_numerics.(name) = vertcat(new_numerics.(name), var(curr_for_dim{:}).(name));
                        end
                        
                        n = length(v_sym);
                        indices = (n_new_sym+1):(n_new_sym+n);
                        curr_for_dim_adj = num2cell(curr(1:dim)+1);
                        var.indices{curr_for_dim_adj{:}} = indices;
                    end
                end
            end
            obj.sym = new_sym;
            obj.numerical_vectors = new_numerics;
        end

        function add_variable_group(obj, name, vars, varargin)
            % Adds a :class:`vdx.VariableGroup` to this vector
            if isfield(obj.variables,name)
                error('Variable or VariableGroup with this name already exists')
            else
                obj.variables.(name) = vdx.VariableGroup(vars, varargin{:});
            end
        end
    end
    
    methods (Access=protected)
        function varargout = dotReference(obj,index_op)
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
                obj.numerical_vectors.(index_op) = varargin{1};
                return
            end
            if isscalar(index_op)
                if ~isfield(obj.variables,name) % Workaround for scalar variables because matlab throws a fit if you try x() = 1;
                    var = vdx.Variable(obj);
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
                var = vdx.Variable(obj);
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
                % Solver needs to be cleared.
                cp.variables.(var_names{ii}).solver = [];
            end
        end
    end

    methods (Access={?vdx.Variable, ?vdx.VariableGroup})
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
                numeric_vals.(name) = p.Results.(name);
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
                    numeric_vals.(name) = numeric_vals.(name) * ones(n,1);
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
                if isa(sym, ['casadi.' obj.casadi_type])
                    % do nothing, directly pass through
                elseif iscell(sym) && length(sym) == 2 &&...
                        ischar(sym{1}) && length(sym{2}) == 1 && isnumeric(sym{2}) && round(sym{2}) == sym{2}
                    sym = define_casadi_symbolic(obj.casadi_type, sym{1}, sym{2});
                else
                    error("Incorrect type")
                end
            else
                if isa(sym, ['casadi.' obj.casadi_type])
                    if sym.is_symbolic % if we can make the names of symbolics nice. However this may be bad if one wants to add single variable constraints in g instead of w.
                        name = split(sym(1).name, '_');
                        name = [name{1:end-1} index_string(p.Results.index)];
                        sym = define_casadi_symbolic(obj.casadi_type, name, size(sym, 1));
                    else
                        % pass through and do nothing
                    end
                elseif iscell(sym) && length(sym) >= 2
                    if ischar(sym{1}) && length(sym{2}) == 1 && round(sym{2}) == sym{2}
                        name = [sym{1} index_string(p.Results.index)];
                        sym = define_casadi_symbolic(obj.casadi_type, name, sym{2});
                    elseif isa(sym{1}, 'casadi.Function')
                        arg_group = vdx.VariableGroup(sym{2}, varargs);
                        fun = sym{1};
                        inds = num2cell(p.Results.index);
                        sym = fun(arg_group(inds{:}));
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
    for i=indices
        istring = [istring '_' num2str(i)];
    end
end
