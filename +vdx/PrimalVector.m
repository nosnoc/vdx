classdef PrimalVector < vdx.Vector
% A subclass of vdx.Vector for primal variables.
%
% :param vdx.Problem problem: Problem which this vector is a member of.
% :param string casadi_type: either 'SX' (default) or 'MX' which determines the kind of CasADi symbolic stored.
% :ivar double lb: Lower bound.
% :ivar double ub: upper bound.
% :ivar double init: Initial Value.
% :ivar double init_mult: Initial value of bound multipliers.
% :ivar double res: Output: primal results.
% :ivar double violation: Output: violation of bounds.
% :ivar double mult: Output: multiplier results for the bounds.
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
