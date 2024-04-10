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
    end

    properties (Access=protected)
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
        allow_non_symbolic_assignment
    end

    methods (Access=public)
        function obj = Vector(problem, varargin)
            p = inputParser;
            addRequired(p, 'problem');
            for name=obj.numerical_properties
                if isfield(obj.numerical_defaults, name)
                    default = obj.numerical_defaults.(name)
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
            

            % Populate defaults
            for name=obj.numerical_properties
                obj.numerical_defaults.(name) = p.Results.(['default_' char(name)]);
            end

            % Populate vectors
            for name=obj.numerical_properties
                obj.numerical_defaults.(name) = [];
            end
            for name=obj.numerical_outputs
                obj.numerical_defaults.(name) = [];
            end
        end

        function print(obj, varargin)
        % Pretty prints this vector with the specified columns.
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
            fprintf(header);

            % iterate over all requested values
            n = size(obj.sym, 1);
            output = [];
            for ii=1:n
                pline = [num2str(ii) '\t\t'];
                if lb
                    pline = [pline sprintf('%-8.5g\t', obj.lb(ii))];
                end
                if ub
                    pline = [pline sprintf('%-8.5g\t', obj.ub(ii))];
                end
                if init
                    pline = [pline sprintf('%-8.5g\t', obj.init(ii))];
                end
                if res
                    pline = [pline sprintf('%-8.5g\t', obj.res(ii))];
                end
                if mult
                    pline = [pline sprintf('%-8.5g\t', obj.mult(ii))];
                end
                if sym
                    pline = [pline char(formattedDisplayText(obj.sym(ii)))];
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
            % TODO(@anton) do we want to also re-organize mult and res?
            new_sym = [];
            new_lb = [];
            new_ub = [];
            new_init = [];

            % First re-normalize 0 dimensional vars (i.e indicies 1x1)
            d_vars = vars_by_depth{1};
            for jj=1:numel(d_vars)
                var = obj.variables.(d_vars{jj});
                n_new_sym = length(new_sym);
                v_sym = var(0);
                v_lb = var(0).lb;
                v_ub = var(0).ub;
                v_init = var(0).init;

                new_sym = vertcat(new_sym, v_sym);
                new_lb = [new_lb; v_lb];
                new_ub = [new_ub; v_ub];
                new_init = [new_init; v_init];

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
                        v_lb = var(curr_for_dim{:}).lb;
                        v_ub = var(curr_for_dim{:}).ub;
                        v_init = var(curr_for_dim{:}).init;

                        new_sym = vertcat(new_sym, v_sym);
                        new_lb = [new_lb; v_lb];
                        new_ub = [new_ub; v_ub];
                        new_init = [new_init; v_init];

                        n = length(v_sym);
                        indices = (n_new_sym+1):(n_new_sym+n);
                        curr_for_dim_adj = num2cell(curr(1:dim)+1);
                        var.indices{curr_for_dim_adj{:}} = indices;
                    end
                end
            end
            obj.sym = new_sym;
            obj.lb = new_lb;
            obj.ub = new_ub;
            obj.init = new_init;
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
        function indices = add_variable(obj, symbolic, varargin)
        % Adds a :class:`vdx.Variable` to the internal symbolic and numeric vectors.
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            for name=obj.numerical_properties
                default = obj.numerical_defaults.(name)
                addParameter(p, name, default);
            end
            parse(p, obj, symbolic, varargin{:});

            symbolic = p.Results.symbolic;

            % Check that symbolic is valid
            if iscell(symbolic) && ~isa(symbolic{1}, 'casadi.Function') && ~obj.allow_nonsymbolic_assignment
                % TODO better error
                error('This vector of class ' class(obj) ' does not allow for {name, size} form of assignment.')
            end
            % TODO more checks here
            
            % Handle non-symbolic input as (name, size) pair
            if iscell(symbolic)
                if ischar(symbolic{1})
                    name = symbolic{1};
                    len = symbolic{2};
                    symbolic = define_casadi_symbolic(class(obj.sym), name, len);
                end
            end

            % Get size and populate possibly empty values
            n = size(symbolic, 1);
            if isempty(lb)
                lb = obj.default_lb*ones(n,1);
            elseif length(lb) == 1
                lb = lb*ones(n,1);
            end
            
            if isempty(ub)
                ub = obj.default_ub*ones(n,1);
            elseif length(ub) == 1
                ub = ub*ones(n,1);
            end

            if isempty(initial)
                initial = obj.default_init*ones(n,1);
            elseif length(initial) == 1
                initial = initial*ones(n,1);
            end
            
            lens = [size(symbolic,1),size(lb,1),size(ub,1), size(initial,1)];
            if ~all(lens == lens(1))
                % TODO(@anton) better error message
                error("mismatched dims")
            end

            n_sym = size(obj.sym, 1);

            obj.sym = vertcat(obj.sym, symbolic);
            obj.lb = [obj.lb; lb];
            obj.ub = [obj.ub; ub];
            obj.init = [obj.init; initial];

            % initialize results and multipliers to zero
            % TODO(@anton) is there a better descision than this?
            obj.res = [obj.res; zeros(n,1)];
            obj.mult = [obj.mult; zeros(n,1)];

            indices = (n_sym+1):(n_sym+n);
        end
    end

    methods (Access=private)
        function sym = eval_symbolic(sym, indices)
        % process sym into a casadi symbolic possibly renaming using indices in the process
        end
    end
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
