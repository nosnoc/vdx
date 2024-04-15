classdef Mpcc < vdx.Problem
    properties (Access=public)
        G vdx.ConstraintVector
        H vdx.ConstraintVector
    end

    methods (Access=public)
        function obj = Mpcc()
            obj = obj@vdx.Problem();
            obj.G = vdx.ConstraintVector(obj);
            obj.H = vdx.ConstraintVector(obj);
        end

        function create_solver(obj, options, plugin)
            w = obj.w(:);
            g = obj.g(:);
            G = obj.G(:);
            H = obj.H(:);
            p = obj.p(:);
            f = obj.f;

            if ~exist('plugin')
                plugin = 'relaxation';
            end

            mpcc_struct = struct('x', w, 'g', g, 'p', p, 'G', G, 'H', H, 'f', f);
            
            % TODO(@anton) figure out if we want to create a standard repository for mpcc solvers or if they should live in nosnoc.
            %              My current thought is to move it out so we don't have circular dependency here. (alternatively move this into nosnoc?)
            obj.solver = nosnoc.solver.mpccsol('name', 'relaxation', mpcc_struct, options);
        end

        function stats = solve(obj)
            mpcc_results = obj.solver('x0', obj.w.init,...
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
            obj.w.res = full(mpcc_results.x);
            obj.w.mult = full(mpcc_results.lam_x);
            obj.g.res = full(mpcc_results.g);
            obj.g.mult = full(mpcc_results.lam_g);
            obj.p.mult = full(mpcc_results.lam_p);
            obj.f_result = full(mpcc_results.f);
            obj.G.res = full(mpcc_results.G);
            obj.H.res = full(mpcc_results.H);
            % TODO(@anton) multipliers for G and H
            stats = obj.solver.stats;
        end

        function mpcc_struct = to_casadi_struct(obj)
            mpcc_struct = struct;
            mpcc_struct.x = mpcc.w.sym;
            mpcc_struct.g = mpcc.g.sym;
            mpcc_struct.p = mpcc.p.sym;
            mpcc_struct.G = mpcc.G.sym;
            mpcc_struct.H = mpcc.H.sym;
            mpcc_struct.f = mpcc.f;
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
