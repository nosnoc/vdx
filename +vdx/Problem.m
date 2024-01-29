classdef Problem < handle &...
        matlab.mixin.Copyable
    properties (Access=public)
        % Primal variabiles
        w
        % Constraints
        g
        % Parameters
        p
        % Objective
        f
    end
    properties (Access=public, NonCopyable)
        % Solver attached to problem
        solver
    end
    properties (Access=private)
        % Variable fun
        variables struct
        constraints struct
    end

    methods (Access=public)
        function obj = Problem()
            obj.w = vdx.Vector(obj, -inf, inf, 0);
            obj.p = vdx.Vector(obj, -inf, inf, 0);
            obj.g = vdx.Vector(obj, 0, 0, 0);
            obj.f = 0;
        end

        function create_solver(obj, casadi_options)
            w = obj.w(:);
            g = obj.g(:);
            p = obj.p(:);
            f = obj.f;

            casadi_nlp = struct('x', w, 'g', g, 'p', p, 'f', f);
            % TODO(@anton) more options than just ipopt.
            % TODO(@anton) options struct should probably be separated into top level and casad opts.
            obj.solver = casadi.nlpsol('proj_fesd', 'ipopt', casadi_nlp, casadi_options);
        end

        function stats = solve(obj)
            nlp_results = obj.solver('x0', obj.w.init,...
                'lbx', obj.w.lb,...
                'ubx', obj.w.ub,...
                'lbg', obj.g.lb,...
                'ubg', obj.g.ub,...
                'lam_g0', obj.g.mult,...% TODO(@anton) perhaps we use init instead of mult.
                'lam_x0', obj.w.mult,...
                'p', obj.p.init);
            if ~obj.solver.stats.success
                %warning("failed to converge")
            end
            obj.w.res = full(nlp_results.x);
            obj.w.mult = full(nlp_results.lam_x);
            obj.g.res = full(nlp_results.g);
            obj.g.mult = full(nlp_results.lam_g);
            obj.p.mult = full(nlp_results.lam_p);
            stats = obj.solver.stats;
        end
    end
    
    methods (Access=protected)
        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);

            cp.w = copy(obj.w);
            cp.w.problem = cp;
            cp.g = copy(obj.g);
            cp.g.problem = cp;
            cp.p = copy(obj.p);
            cp.p.problem = cp;
            cp.f = obj.f;
        end
    end
end
