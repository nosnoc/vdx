classdef ParameterVector < vdx.Vector
    properties (Access=protected)
        numerical_defaults = struct('lb', -inf, 'ub', inf, 'init', 0);
    end
    properties (Constant, Hidden)
        numerical_properties = ["val"];
        numerical_outputs = ["lambda"];

        allow_nonscalar_symbolics = false;
        allow_nonsymbolic_assignment = true;
        % TODO also do bound violation in output
    end

    methods
        function obj = ParameterVector(problem, varargin)
            obj = obj@vdx.Vector(problem, varargin{:});
        end
    end
end
