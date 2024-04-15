classdef PrimalVector < vdx.Vector
    properties (Access=protected)
        numerical_defaults = struct('lb', -inf, 'ub', inf, 'init', 0);
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
    end
end
