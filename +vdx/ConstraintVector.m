classdef ConstraintVector < vdx.Vector
    properties (Access=protected)
        numerical_defaults = struct('lb', 0, 'ub', 0, 'init_lambda', 0);
    end
    properties (Constant, Hidden)
        numerical_properties = ["lb", "ub", "init_lambda"];
        numerical_outputs = ["eval", "violation", "lambda"];

        allow_nonscalar_symbolics = true;
        allow_nonsymbolic_assignment = false;
        % TODO also do bound violation in output
    end

    methods
        function obj = PrimalVector(problem, varargin)
            obj = obj@vdx.Vector(problem, varargin{:});
        end
    end
end
