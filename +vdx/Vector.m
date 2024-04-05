classdef Vector < handle &...
        matlab.mixin.indexing.RedefinesDot &...
        matlab.mixin.indexing.RedefinesParen &...
        dynamicprops &...
        matlab.mixin.Copyable
    properties (Access=public)
        % Symbolic vector that this wraps TODO(@anton) possibly rename?
        sym
        % Numeric vectors containing lower and upper bounds and initial data
        lb
        ub
        init
        % Numeric vectors for results
        res
        mult
        % Default values for bounds and init
        default_lb
        default_ub
        default_init
    end
    properties (Access=private)
        % pointer to parent problem
        % TODO(@anton) think about this and possibly provide functionality to move between problems. (low priority)
        problem
        % casadi type
        %casadi_type
        % Internal struct of index tracking variables
        variables struct
    end

    properties (Dependent)
        
    end

    methods (Access=public)
        function obj = Vector(problem, varargin)
            p = inputParser;
            addRequired(p, 'problem');
            addOptional(p, 'lb', -inf);
            addOptional(p, 'ub', inf);
            addOptional(p, 'init', 0);
            addParameter(p, 'casadi_type', 'SX');
            parse(p, problem, varargin{:});
            
            obj.problem = problem;
            obj.variables = struct;

            % Populate defaults
            obj.default_lb = p.Results.lb;
            obj.default_ub = p.Results.ub;
            obj.default_init = p.Results.init;

            sym = casadi.(p.Results.casadi_type);
            obj.lb = [];
            obj.ub = [];
            obj.init = [];
        end
        
        function indices = add_variable(obj, symbolic, varargin)
        % ADD_VARIABLE  Adds a variable to the internal symbolic and numeric vectors.
        %   TODO(@anton) Perhaps this should be private and renamed. (medium priority)
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            addOptional(p, 'lb', []);
            addOptional(p, 'ub', []);
            addOptional(p, 'initial', []);
            parse(p, obj, symbolic, varargin{:});

            symbolic = p.Results.symbolic;
            lb = p.Results.lb;
            ub = p.Results.ub;
            initial = p.Results.initial;

            % Handle non-symbolic input as (name, size) pair
            if iscell(symbolic)
                if ischar(symbolic{1})
                    name = symbolic{1};
                    len = symbolic{2};
                    symbolic = define_casadi_symbolic(class(sym), name, len);
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

        function add_variable_group(obj, name, vars, varargin)
            if isfield(obj.variables,name)
                error('Variable or VariableGroup with this name already exists')
            else
                obj.variables.(name) = vdx.VariableGroup(vars, varargin{:});
            end
        end

        function out = cat(dim,varargin)
            error('Concatenation not (yet) supported')
            % TODO(@anton) This is certainly possible but will take some work 
        end

        function varargout = size(obj,varargin)
            varargout = size(sym, varargin{:});
            %TODO(anton) needs to return correct values for varargin and perhaps other cases?
        end

        function print(obj, varargin)
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
            n = size(sym, 1);
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
                    pline = [pline char(formattedDisplayText(sym(ii)))];
                end
                pline = [pline, '\n'];
                output = [output pline];
            end
            fprintf(output);
        end

        function sort_by_index(obj)
        % SORT_BY_INDEX Sorts this vector so that the vectors occur in column major order with lower dimensional variables
        %               always occuring before higher dimensional variables. This can be useful to recover any structure in the 
        %               problem that comes from the structure of constaraints and variables.
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
                %disp(dims_to_process)
                %disp(curr)
                for dim=dims_to_process
                    d_vars = vars_by_depth{dim+1};
                    curr_for_dim = num2cell(curr(1:dim));
                    %disp(curr(1:dim))
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
            sym = new_sym;
            obj.lb = new_lb;
            obj.ub = new_ub;
            obj.init = new_init;
        end
    end
    
    methods (Access=protected)
        % Dot reference overrides
        function varargout = dotReference(obj,index_op)
            name = index_op(1).Name;
            if ~isfield(obj.variables, name)
                var = vdx.Variable(obj,[]);
                obj.variables.(name) = var; % TODO(@anton) rename this. Talk to Armin.
                P = obj.addprop(name);
                obj.(name) = var;
            end
            varargout{1} = obj.variables.(index_op);
        end

        function obj = dotAssign(obj,index_op,varargin)
            if string(version('-release')) < "2023a"
                persistent matlab_bug_workaround_to_prevent_garbage_collection;
                matlab_bug_workaround_to_prevent_garbage_collection = varargin{1};
            end
            name = index_op(1).Name;
            if ~isfield(obj.variables,name)
                var = vdx.Variable(obj,[]);
                obj.variables.(name) = var; % TODO(@anton) rename this. Talk to Armin.
                P = obj.addprop(name);
                obj.(name) = var;
            end
            obj.variables.(index_op) = varargin{1};
        end
        
        function n = dotListLength(obj,index_op,indexContext)
            n=1;
        end

        function varargout = parenReference(obj, index_op)
            varargout{1} = sym.(index_op); % TODO(@anton) this should be sufficient to pass through
        end

        function obj = parenAssign(obj,index_op,varargin)
            error("Raw parenAssign is likely an error.")
            % TODO(@anton) is it?
        end

        function obj = parenDelete(obj,index_op)
            error('Raw parenDelete is unsupported')
            % TODO(@anton) is it?
        end

        function n = parenListLength(obj,index_op,ctx)
           n = 1;
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
