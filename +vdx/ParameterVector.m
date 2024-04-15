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
    end
end
