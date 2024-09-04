classdef ParameterVector < vdx.Vector
% A subclass of vdx.Vector for parameters.
%
% :param vdx.Problem problem: Problem which this vector is a member of.
% :param string casadi_type: either 'SX' (default) or 'MX' which determines the kind of CasADi symbolic stored.
% :ivar double val: Value of the parameter.
% :ivar double mult: Output: multiplier (sensitivity) of the problem to the parameter value. 
    properties (Access=protected)
        numerical_defaults = struct('val', 0);
    end
    properties (Constant, Hidden)
        numerical_properties = ["val"];
        numerical_outputs = ["mult"];

        allow_nonscalar_symbolics = false;
        allow_nonsymbolic_assignment = true;
        % TODO also do bound violation in output
    end

    methods
        function obj = ParameterVector(problem, varargin)
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
            vec_struct.numerical_vectors.val = obj.numerical_vectors.val;
            vec_struct.numerical_vectors.mult = obj.numerical_vectors.mult;
            
            json = jsonencode(vec_struct, varargin{:});
        end
    end

    methods(Static)
        function vec = from_json(vec_struct)
            vec = vdx.ParameterVector([]);
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
