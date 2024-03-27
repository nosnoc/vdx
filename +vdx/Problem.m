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
        % objective value
        f_result
    end
    properties (Access=public, NonCopyable)
        % Solver attached to problem
        solver
    end

    methods (Access=public)
        function obj = Problem(varargin)
            p = inputParser;
            addParameter(p, 'casadi_type', 'SX');
            parse(p, varargin{:});
            
            obj.w = vdx.Vector(obj, -inf, inf, 0, 'casadi_type', p.Results.casadi_type);
            obj.p = vdx.Vector(obj, -inf, inf, 0, 'casadi_type', p.Results.casadi_type);
            obj.g = vdx.Vector(obj, 0, 0, 0, 'casadi_type', p.Results.casadi_type);
            obj.f = 0;
            obj.f_result = 0;
        end

        function create_solver(obj, casadi_options, plugin)
            w = obj.w(:);
            g = obj.g(:);
            p = obj.p(:);
            f = obj.f;

            if ~exist('plugin')
                plugin = 'ipopt';
            end

            casadi_nlp = struct('x', w, 'g', g, 'p', p, 'f', f);
            obj.solver = casadi.nlpsol('proj_fesd', plugin, casadi_nlp, casadi_options);
        end

        function [stats, nlp_results] = solve(obj)
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
            obj.f_result = full(nlp_results.f);
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
