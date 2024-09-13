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
    end

    methods
        function obj = ConstraintVector(problem, varargin)
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

            fun = casadi.Function('foo', {obj.problem.w.sym, obj.problem.p.sym}, {obj.sym});
            
            vec_struct.sym = fun.serialize;
            vec_struct.casadi_type = obj.casadi_type;
            vec_struct.len = obj.len;
            vec_struct.numerical_defaults = obj.numerical_defaults;
            vec_struct.numerical_vectors.lb = obj.numerical_vectors.lb;
            vec_struct.numerical_vectors.ub = obj.numerical_vectors.ub;
            vec_struct.numerical_vectors.init_mult = obj.numerical_vectors.init_mult;
            vec_struct.numerical_vectors.eval = obj.numerical_vectors.eval;
            vec_struct.numerical_vectors.violation = obj.numerical_vectors.violation;
            vec_struct.numerical_vectors.mult = obj.numerical_vectors.mult;
            
            json = jsonencode(vec_struct, varargin{:});
        end
    end

    methods(Static)
        function vec = from_json(vec_struct, problem)
            vec = vdx.ConstraintVector([]);
            vec.casadi_type = vec_struct.casadi_type;
            vec.len = vec_struct.len;
            vec.numerical_defaults = vec_struct.numerical_defaults;
            vec.numerical_vectors = vec_struct.numerical_vectors;

            fun = casadi.Function.deserialize(vec_struct.sym);
            vec.sym = fun(problem.w.sym, problem.p.sym);

            names = fields(vec_struct.variables);
            for ii=1:numel(names)
                name = names{ii};
                vec.variables.(name) = vdx.Variable.from_json(vec_struct.variables.(name));
                vec.variables.(name).vector = vec;
                vec.addprop(name);
                vec.(name) = vec.variables.(name);
            end
        end
    end
    
    methods (Access={?vdx.Variable, ?vdx.VariableGroup, ?vdx.ConstraintVector})
        function indices = add_variable(obj, indices, symbolic, varargin)
            import casadi.*
            % Just add variable if no relaxation passed
            if ~length(varargin) || class(varargin{end}) ~= "vdx.RelaxationStruct"
                indices = add_variable@vdx.Vector(obj, indices, symbolic, varargin{:});
                return
            end
            relax = varargin{end};
            % Exit early if relaxation mode none
            if relax.mode == vdx.RelaxationMode.NONE
                indices = add_variable@vdx.Vector(obj, indices, symbolic, varargin{1:end-1});
                return
            end

            % Build and parse inputs
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

            % get weight parameter
            if ~obj.problem.p.has_var(relax.weight_name);
                obj.problem.p.(relax.weight_name) = {{relax.weight_name, 1}, 1};
            end
            weight = obj.problem.p.(relax.weight_name)();

            lb = numeric_vals.lb;
            ub = numeric_vals.ub;
            new_sym = {};
            new_lb = [];
            new_ub = [];
            % TODO not just SX and names
            slacks = {};
            switch relax.mode
              case vdx.RelaxationMode.ELL_1 % l_1
                for ii=1:n
                    if lb(ii) == ub(ii) % equality
                        slack = SX.sym(relax.slack_name, 2);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)-lb(ii)-slack(1) + slack(2)];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;0];
                    elseif lb(ii) >-inf & ub(ii) < inf % twosided
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)-lb(ii)+slack;
                            -(symbolic(ii)-ub(ii))+slack];
                        new_lb = [new_lb;0;0];
                        new_ub = [new_ub;inf;inf];
                    elseif lb(ii) == -inf % only ub
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [-(symbolic(ii)-ub(ii))+slack];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    else % only lb
                        slack = SX.sym(relax.slack_name, 1);
                        slacks{ii} = slack;
                        new_sym{ii} = [symbolic(ii)-lb(ii)+slack];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    end
                end
                slacks = vertcat(slacks{:});
                obj.problem.w.(relax.slack_name)(indices{:}) = {slacks, 0, inf, 100};
                new_sym = vertcat(new_sym{:});
                obj.problem.f = obj.problem.f + weight*sum(slacks); % TODO weight

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
                obj.problem.f = obj.problem.f + weight*sum(slacks.^2); % TODO weight

                indices = add_variable@vdx.Vector(obj, indices, new_sym, numeric_vals.lb, numeric_vals.ub, numeric_vals.init_mult);
              case vdx.RelaxationMode.ELL_INF % l_inf
                if ~obj.problem.w.has_var(relax.slack_name);
                    % create slacks in the form (s_eq, s^+, s^-, s_max)
                    obj.problem.w.(relax.slack_name) = {{relax.slack_name, 4}, [-inf;0;0;0], [inf;inf;inf;inf], [0; 100; 100; 100]};
                    slack = obj.problem.w.(relax.slack_name)();
                    % Add abs constraint s_eq = s^+ - s^-, s_max >= s^+, s_max >= s^-, and penalty on s_max
                    obj.problem.g.(['ell_inf_' relax.slack_name]) = {[slack(1) - slack(2) + slack(3);
                                                                       slack(4) - slack(2);
                                                                       slack(4) - slack(3)], [0;0;0], [0;inf;inf]};
                    obj.problem.f = obj.problem.f + weight*slack(4);
                else
                    slack = obj.problem.w.(relax.slack_name)();
                end

                s_eq = slack(1);
                s = slack(4);

                for ii=1:n
                    if lb(ii) == ub(ii) % equality
                        new_sym{ii} = [symbolic(ii)-lb(ii)+s_eq];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;0];
                    elseif lb(ii) >-inf & ub(ii) < inf % twosided
                        new_sym{ii} = [symbolic(ii)-lb(ii)+s;
                            -(symbolic(ii)-ub(ii))+s];
                        new_lb = [new_lb;0;0];
                        new_ub = [new_ub;inf;inf];
                    elseif lb(ii) == -inf % only ub
                        new_sym{ii} = [-(symbolic(ii)-ub(ii))+s];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    else % only lb
                        new_sym{ii} = [symbolic(ii)-lb(ii)+s];
                        new_lb = [new_lb;0];
                        new_ub = [new_ub;inf];
                    end
                end

                new_sym = vertcat(new_sym{:});
                indices = add_variable@vdx.Vector(obj, indices, new_sym, new_lb, new_ub, numeric_vals.init_mult);
            end
        end
    end
end
