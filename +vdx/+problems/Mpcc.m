classdef Mpcc < vdx.Problem
    properties (Access=public)
        G vdx.ConstraintVector
        H vdx.ConstraintVector

        mpcc_results
    end

    methods (Access=public)
        function obj = Mpcc()
            obj = obj@vdx.Problem();
            obj.G = vdx.ConstraintVector(obj);
            obj.H = vdx.ConstraintVector(obj);
        end

        function create_solver(obj, options, plugin)
            if ~exist('plugin')
                plugin = 'scholtes_ineq';
            end

            %mpcc_struct = obj.to_casadi_struct();
            
            % TODO(@anton) figure out if we want to create a standard repository for mpcc solvers or if they should live in nosnoc.
            %              My current thought is to move it out so we don't have circular dependency here. (alternatively move this into nosnoc?)
            obj.solver = nosnoc.solver.mpccsol('Mpcc solver', plugin, obj.to_casadi_struct(), options);
        end

        function stats = solve(obj)
            mpcc_results = obj.solver('x0', obj.w.init,...
                'lbx', obj.w.lb,...
                'ubx', obj.w.ub,...
                'lbg', obj.g.lb,...
                'ubg', obj.g.ub,...
                'lam_g0', obj.g.init_mult,...% TODO(@anton) perhaps we use init instead of mult.
                'lam_x0', obj.w.init_mult,...
                'p', obj.p.val);
            if ~obj.solver.stats.success
                %warning("failed to converge")
            end
            obj.w.res = full(mpcc_results.x);
            obj.w.mult = full(mpcc_results.lam_x);
            obj.g.eval = full(mpcc_results.g);
            obj.g.mult = full(mpcc_results.lam_g);
            % TODO(@anton) can we figure out parameter multipliers correctly from homotopy?
            %obj.p.mult = full(mpcc_results.lam_p);
            obj.f_result = full(mpcc_results.f);
            obj.G.eval = full(mpcc_results.G);
            obj.H.eval = full(mpcc_results.H);
            % TODO(@anton) multipliers for G and H

            % Calculate violations:
            w_lb_viol = max(obj.w.lb - obj.w.res, 0);
            w_ub_viol = max(obj.w.res - obj.w.ub, 0);
            obj.w.violation = max(w_lb_viol, w_ub_viol);
            g_lb_viol = max(obj.g.lb - obj.g.eval, 0);
            g_ub_viol = max(obj.g.eval - obj.g.ub, 0);
            obj.g.violation = max(g_lb_viol, g_ub_viol);

            stats = obj.solver.stats;
            obj.mpcc_results = mpcc_results;
        end

        function mpcc_struct = to_casadi_struct(obj)
            mpcc_struct = struct;
            mpcc_struct.x = obj.w.sym;
            mpcc_struct.g = obj.g.sym;
            mpcc_struct.p = obj.p.sym;
            mpcc_struct.G = obj.G.sym;
            mpcc_struct.H = obj.H.sym;
            mpcc_struct.f = obj.f;
        end

        function print(obj, varargin)
            w_out = obj.w.to_string(varargin{:});
            p_out = obj.p.to_string(varargin{:});
            g_out = obj.g.to_string(varargin{:});
            G_out = obj.G.to_string(varargin{:});
            H_out = obj.H.to_string(varargin{:});

            n_longest = max([longest_line(w_out),
                              longest_line(p_out),
                              longest_line(g_out),
                              longest_line(G_out),
                              longest_line(H_out)]);

            hline(1:n_longest) = '-';
            hline = [hline '\n'];
            
            fprintf(hline);
            fprintf("Primal Variables\n");
            fprintf(hline);
            fprintf(w_out);
            fprintf(hline);
            fprintf("Parameters\n");
            fprintf(hline);
            fprintf(p_out);
            fprintf(hline);
            fprintf("Constraints\n");
            fprintf(hline);
            fprintf(g_out);
            fprintf(hline);
            fprintf("G Complementarity\n");
            fprintf(hline);
            fprintf(G_out);
            fprintf(hline);
            fprintf("H Complementarity\n");
            fprintf(hline);
            fprintf(H_out);
            fprintf(hline);
            fprintf("Objective\n");
            fprintf(hline);
            print_casadi_vector(obj.f);
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
