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

        function init_struct = to_solver_initialization(obj)
            init_struct = struct;
            init_struct.x0 = obj.w.init;
            init_struct.lbx = obj.w.lb;
            init_struct.ubx = obj.w.ub;
            init_struct.lam_x0 = obj.w.init_mult;
            init_struct.lbg = obj.g.lb;
            init_struct.ubg = obj.g.ub;
            init_struct.lam_g0 = obj.g.init_mult;
            init_struct.p0 = obj.p.val;
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

        function json = jsonencode(obj, varargin)
            obj.finalize_assignments();
            idx_struct = struct();

            idx_struct.w = obj.w;
            idx_struct.g = obj.g;
            idx_struct.p = obj.p;
            idx_struct.G = obj.G;
            idx_struct.H = obj.H;

            f_fun = casadi.Function('f', {obj.w.sym, obj.p.sym}, {obj.f});
            idx_struct.f = f_fun.serialize;
            
            json = jsonencode(idx_struct, varargin{:});
        end

        function json = to_casadi_json(obj)
            obj.finalize_assignments();
            mpcc_struct = obj.to_casadi_struct();

            json_struct.w = mpcc_struct.x.serialize();
            json_struct.lbw = obj.w.lb;
            json_struct.ubw = obj.w.ub;
            json_struct.w0 = obj.w.init;
            json_struct.p = mpcc_struct.p.serialize();
            json_struct.p0 = obj.p.val;
            json_struct.g_fun = casadi.Function('g', {mpcc_struct.x,mpcc_struct.p}, {mpcc_struct.g}).serialize();
            json_struct.lbw = obj.g.lb;
            json_struct.ubw = obj.g.ub;
            json_struct.G_fun = casadi.Function('G', {mpcc_struct.x,mpcc_struct.p}, {mpcc_struct.G}).serialize();
            json_struct.H_fun = casadi.Function('H', {mpcc_struct.x,mpcc_struct.p}, {mpcc_struct.H}).serialize();
            json_struct.f = casadi.Function('f', {mpcc_struct.x,mpcc_struct.p}, {mpcc_struct.f}).serialize();

            json = jsonencode(json_struct, "ConvertInfAndNaN", false, "PrettyPrint", true);
        end

        function finalize_assignments(obj)
            obj.w.apply_queued_assignments;
            obj.g.apply_queued_assignments;
            obj.p.apply_queued_assignments;
            obj.G.apply_queued_assignments;
            obj.H.apply_queued_assignments;
        end
    end

    methods(Static)
        function mpcc = from_json(txt)
            mpcc_struct = jsondecode(txt);
            mpcc = vdx.problems.Mpcc();
            mpcc.w = vdx.PrimalVector.from_json(mpcc_struct.w);
            mpcc.w.problem = mpcc;
            mpcc.p = vdx.ParameterVector.from_json(mpcc_struct.p);
            mpcc.p.problem = mpcc;
            mpcc.g = vdx.ConstraintVector.from_json(mpcc_struct.g, mpcc);
            mpcc.g.problem = mpcc;
            mpcc.G = vdx.ConstraintVector.from_json(mpcc_struct.G, mpcc);
            mpcc.G.problem = mpcc;
            mpcc.H = vdx.ConstraintVector.from_json(mpcc_struct.H, mpcc);
            mpcc.H.problem = mpcc;

            f_fun = casadi.Function.deserialize(mpcc_struct.f);
            mpcc.f = f_fun(mpcc.w.sym, mpcc.p.sym);
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
