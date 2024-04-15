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
    end
end
