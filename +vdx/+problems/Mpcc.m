classdef Mpcc < vdx.Problem
    properties (Access=public)
        G
        H
    end

    methods (Access=public)
        function obj = Mpcc()
            obj = obj@vdx.Problem();
            obj.G = vdx.Vector(obj, 0, inf, 0);
            obj.H = vdx.Vector(obj, 0, inf, 0);
        end

        function create_solver(obj, options)
            w = obj.w(:);
            g = obj.g(:);
            G = obj.G(:);
            H = obj.H(:);
            p = obj.p(:);
            f = obj.f;

            mpcc_struct = struct('x', w, 'g', g, 'p', p, 'G', G, 'H', H, 'f', f);
            
            % TODO(@anton) figure out if we want to create a standard repository for mpcc solvers or if they should live in nosnoc
            obj.solver = nosnoc.solver.mpccsol('name', 'type', mpcc_struct, options); 
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
            cp.G = copy(obj.G);
            cp.G.problem = cp;
            cp.H = copy(obj.H);
            cp.H.problem = cp;
            cp.p = copy(obj.p);
            cp.p.problem = cp;
            cp.f = obj.f;
        end
    end
end
