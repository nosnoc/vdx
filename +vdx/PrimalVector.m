classdef PrimalVector < vdx.Vector
% A subclass of vdx.Vector for primal variables.
%
%:param vdx.Problem problem: Problem which this vector is a member of.
%:param string casadi_type: either 'SX' (default) or 'MX' which determines the kind of CasADi symbolic stored.
%
%Numerical Properties:
%    lb (double): Lower bound
%    ub (double): Upper bound
    properties (Access=protected)
        numerical_defaults = struct('lb', -inf, 'ub', inf, 'init', 0, 'init_mult', 0);
    end
    properties (Constant, Hidden)
        numerical_properties = ["lb", "ub", "init", "init_mult"];
        numerical_outputs = ["res", "violation", "mult"];

        allow_nonscalar_symbolics = false;
        allow_nonsymbolic_assignment = true;
        % TODO also do bound violation in output
    end

    methods
        function obj = PrimalVector(problem, varargin)
            obj = obj@vdx.Vector(problem, varargin{:});
        end

        function json = jsonencode(obj, varargin)
            obj.apply_queued_assignments();
            vec_struct = struct();

            names = fields(obj.variables);
            for ii=1:numel(names)
                name = names{ii};
                vec_struct.variables.(name) = obj.variables.(name);
            end

            vec_struct.sym = obj.sym.serialize;
            vec_struct.casadi_type = obj.casadi_type;
            vec_struct.len = obj.len;
            vec_struct.numerical_defaults = obj.numerical_defaults;
            vec_struct.numerical_vectors.lb = obj.numerical_vectors.lb;
            vec_struct.numerical_vectors.ub = obj.numerical_vectors.ub;
            vec_struct.numerical_vectors.init = obj.numerical_vectors.init;
            vec_struct.numerical_vectors.init_mult = obj.numerical_vectors.init_mult;
            vec_struct.numerical_vectors.res = obj.numerical_vectors.res;
            vec_struct.numerical_vectors.violation = obj.numerical_vectors.violation;
            vec_struct.numerical_vectors.mult = obj.numerical_vectors.mult;
            
            json = jsonencode(vec_struct, varargin{:});
        end
    end

    methods(Static)
        function vec = from_json(vec_struct)
            vec = vdx.PrimalVector([]);
            vec.casadi_type = vec_struct.casadi_type;
            vec.len = vec_struct.len;
            vec.numerical_defaults = vec_struct.numerical_defaults;
            vec.numerical_vectors = vec_struct.numerical_vectors;
            vec.sym = casadi.(vec.casadi_type).deserialize(vec_struct.sym);

            names = fields(vec_struct.variables);
            for ii=1:numel(names)
                name = names{ii};
                obj.variables.(name) = vdx.Variable.from_json(vec_struct.variables.(name));
                obj.variables.(name).vector = vec;
            end
        end
    end
end
