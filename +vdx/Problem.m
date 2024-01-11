classdef Problem < handle 
    properties (Access=public)
        % Primal variabiles
        w
        % Constraints
        g
        % Parameters
        p
        % Objective
        f
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
            % TODO(@anton) more options than just ipopt
            obj.solver = casadi.nlpsol('proj_fesd', 'ipopt', casadi_nlp, casadi_options);
        end

        function solve(obj)
            nlp_results = obj.solver('x0', obj.w.init,...
                'lbx', obj.w.lb,...
                'ubx', obj.w.ub,...
                'lbg', obj.g.lb,...
                'ubg', obj.g.ub,...
                'p', obj.p.init);

            obj.w.res = full(nlp_results.x);
            obj.w.mult = full(nlp_results.lam_x);
            obj.g.res = full(nlp_results.g);
            obj.g.mult = full(nlp_results.lam_g);
            obj.p.mult = full(nlp_results.lam_p);
        end
    end
    
    methods (Access=protected)
    end
end
