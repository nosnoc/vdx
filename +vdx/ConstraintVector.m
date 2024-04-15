classdef ConstraintVector < vdx.Vector
% A subclass of vdx.Vector for algebraic constraints.
%
% :param vdx.Problem problem: Problem which this vector is a member of.
% :param string casadi_type: either 'SX' (default) or 'MX' which determines the kind of CasADi symbolic stored.
% :ivar double lb: Lower bound.
% :ivar double ub: Upper bound.
% :ivar double init_mult: Initial value of bound multipliers.
% :ivar double eval: Output: evaluation of constraint functions.
% :ivar double violation: Output: violation of bounds.
% :ivar double mult: Output: multiplier results for the bounds.
    properties (Access=protected)
        numerical_defaults = struct('lb', 0, 'ub', 0, 'init_mult', 0);
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
