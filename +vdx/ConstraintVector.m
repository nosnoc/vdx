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

    methods (Access={?vdx.Variable, ?vdx.VariableGroup, ?vdx.ConstraintVector})
        function indices = add_variable(obj, indices, symbolic, varargin)
            import casadi.*
            if ~length(varargin) || class(varargin{end}) ~= "vdx.RelaxationStruct"
                indices = add_variable@vdx.Vector(obj, indices, symbolic, varargin{:});
                return
            end

            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'symbolic');
            for name=obj.numerical_properties
                default = obj.numerical_defaults.(name);
                addOptional(p, name, default);
            end
            parse(p, obj, symbolic, varargin{1:(end-1)});

            symbolic = p.Results.symbolic;
            numeric_vals = struct;
            for name=obj.numerical_properties
                numeric_vals.(name) = p.Results.(name);
            end
            
            n = size(symbolic, 1);
            for name=obj.numerical_properties
                if length(p.Results.(name)) == 1
                    numeric_vals.(name) = numeric_vals.(name) * ones(n,1);
                end
            end

            lb = numeric_vals.lb;
            ub = numeric_vals.ub;
            relax = varargin{end};
            new_sym = {};
            new_lb = [];
            new_ub = [];
            % TODO not just SX and names
            slacks = {};
            switch relax.mode
              case vdx.RelaxationMode.ELL_1 % l_1
                for ii=1:n
                    if lb(ii) >-inf & ub(ii) < inf % twosided
                        slack = SX.sym(relax.slack_name, 2);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)-lb(ii)-slack(1);
                            -(symbolic(ii)-ub(ii)-slack(2))];
                        new_lb = [new_lb;0;0];
                        new_ub = [new_ub;inf;inf];
                    elseif lb(ii) == -inf % only ub
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [-(symbolic(ii)-ub(ii)-slack)];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    else % only lb
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)-lb(ii)-slack];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    end
                end
                slacks = vertcat(slacks{:});
                obj.problem.w.(relax.slack_name)(indices{:}) = {slacks, 0, inf, 100};
                new_sym = vertcat(new_sym{:});
                obj.problem.f = obj.problem.f + sum(slacks); % TODO weight

                indices = add_variable@vdx.Vector(obj, indices, new_sym, new_lb, new_ub, numeric_vals.init_mult);
              case vdx.RelaxationMode.ELL_2 % l_2
                for ii=1:n
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)+slack];
                end
                slacks = vertcat(slacks{:});
                obj.problem.w.(relax.slack_name)(indices{:}) = {slacks, -inf, inf, 0};
                new_sym = vertcat(new_sym{:});
                obj.problem.f = obj.problem.f + sum(slacks.^2); % TODO weight

                indices = add_variable@vdx.Vector(obj, indices, new_sym, numeric_vals.lb, numeric_vals.ub, numeric_vals.init_mult);
              case vdx.RelaxationMode.ELL_INF % l_inf
                % TODO
            end
        end
    end
end
