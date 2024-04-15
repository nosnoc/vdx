classdef ConstraintVector < vdx.Vector
    properties (Access=protected)
        numerical_defaults = struct('lb', 0, 'ub', 0, 'init_lambda', 0);
    end
    properties (Constant, Hidden)
        numerical_properties = ["lb", "ub", "init_mult"];
        numerical_outputs = ["eval", "violation", "mult"];

        allow_nonscalar_symbolics = true;
        allow_nonsymbolic_assignment = false;
        % TODO also do bound violation in output
    end

    methods
        function obj = ConstraintVector(problem, varargin)
            obj = obj@vdx.Vector(problem, varargin{:});
        end
    end
end
