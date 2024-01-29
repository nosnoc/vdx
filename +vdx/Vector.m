classdef Vector < handle &...
        matlab.mixin.indexing.RedefinesDot &...
        matlab.mixin.indexing.RedefinesParen &...
        matlab.mixin.Copyable
    properties (Access=public)
        % Symbolic vector that this wraps TODO(@anton) possibly rename?
        w
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
        % Internal struct of index tracking variables
        variables struct
    end

    methods (Access=public)
        function obj = Vector(problem, varargin)
            p = inputParser;
            addRequired(p, 'problem');
            addOptional(p, 'lb', -inf);
            addOptional(p, 'ub', inf);
            addOptional(p, 'init', 0);
            parse(p, problem, varargin{:});
            
            obj.problem = problem;
            obj.variables = struct;

            % Populate defaults
            obj.default_lb = p.Results.lb;
            obj.default_ub = p.Results.ub;
            obj.default_init = p.Results.init;

            % TODO(@anton) How we handle MX and SX should possibly be decided here in a smarter way,
            %              maybe as a flag to the constructor
            obj.w = casadi.SX();
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
                name = symbolic{1};
                len = symbolic{2};

                symbolic = define_casadi_symbolic(class(obj.w), name, len);
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

            n_w = size(obj.w, 1);

            obj.w = vertcat(obj.w, symbolic);
            obj.lb = [obj.lb; lb];
            obj.ub = [obj.ub; ub];
            obj.init = [obj.init; initial];

            % initialize results and multipliers to zero
            % TODO(@anton) is there a better descision than this?
            obj.res = [obj.res; zeros(n,1)];
            obj.mult = [obj.mult; zeros(n,1)];

            indices = (n_w+1):(n_w+n);
        end

        function out = cat(dim,varargin)
            error('Concatenation not (yet) supported')
            % TODO(@anton) This is certainly possible but will take some work 
        end

        function varargout = size(obj,varargin)
            varargout = size(obj.w, varargin{:});
            %TODO(anton) needs to return correct values for varargin and perhaps other cases?
        end
    end
    
    methods (Access=protected)
        % Dot reference overrides
        function varargout = dotReference(obj,index_op)
            name = index_op(1).Name;
            if ~isfield(obj.variables, name)
                obj.variables.(name) = vdx.Variable(obj,[]); % TODO(@anton) rename this. Talk to Armin.
            end
            varargout{1} = obj.variables.(index_op);
        end

        function obj = dotAssign(obj,index_op,varargin)
            name = index_op(1).Name;
            if ~isfield(obj.variables,name)
                obj.variables.(name) = vdx.Variable(obj,[]);
            end
            obj.variables.(index_op) = varargin{1};
        end
        
        function n = dotListLength(obj,index_op,indexContext)
            n=1;
        end

        function varargout = parenReference(obj, index_op)
            varargout{1} = obj.w.(index_op); % TODO(@anton) this should be sufficient to pass through
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
