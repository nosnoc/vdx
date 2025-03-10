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
                plugin = 'reg_homotopy';
            end

            %mpcc_struct = obj.to_casadi_struct();
            
            % TODO(@anton) figure out if we want to create a standard repository for mpcc solvers or if they should live in nosnoc.
            %              My current thought is to move it out so we don't have circular dependency here. (alternatively move this into nosnoc?)
            obj.solver = nosnoc.solver.mpccsol('Mpcc solver', plugin, obj, options);
        end

        function stats = solve(obj,params)
            arguments
                obj
                params.IG = []
                params.IH = []
            end
            y0 = params.IH;
            mpcc_results = obj.solver('x0', obj.w.init,...
                'lbx', obj.w.lb,...
                'ubx', obj.w.ub,...
                'lbg', obj.g.lb,...
                'ubg', obj.g.ub,...
                'lam_g0', obj.g.init_mult,...% TODO(@anton) perhaps we use init instead of mult.
                'lam_x0', obj.w.init_mult,...
                'p', obj.p.val,...
                'y0', y0);% TODO(@anton) I actually hate this y0 interface :(
            if ~obj.solver.stats.success
                %warning("failed to converge")
            end
            obj.w.res = full(mpcc_results.x);
            try
                obj.w.mult = full(mpcc_results.lam_x);
            catch
                obj.w.mult = nan(size(obj.w.res));
            end
            try
                obj.g.eval = full(mpcc_results.g);
            catch
                obj.g.eval = nan(size(obj.g.sym));
            end
            try
                obj.g.mult = full(mpcc_results.lam_g);
            catch
                obj.g.mult = nan(size(obj.g.sym));
            end
            % TODO(@anton) can we figure out parameter multipliers correctly from homotopy?
            %obj.p.mult = full(mpcc_results.lam_p)
            obj.f_result = full(mpcc_results.f);
            try
                obj.G.eval = full(mpcc_results.G);
            catch
                obj.G.eval = nan(size(obj.G.sym));
            end
            try
                obj.H.eval = full(mpcc_results.H);
            catch
                obj.H.eval = nan(size(obj.H.sym));
            end
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
            json_struct.lbg = obj.g.lb;
            json_struct.ubg = obj.g.ub;
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

        function [bnlp, solver_initialization] = generate_bnlp_for_active_set(obj, I1, I2)
            import casadi.*
            % TODO: might be worth using vdx to maintain order for the bnlp but that may be too expensive
            % bnlp = vdx.Problem();
            
            % build nlp:
            w = obj.w.sym; lbw = obj.w.lb; ubw = obj.w.ub;
            g = obj.g.sym; lbg = obj.g.lb; ubg = obj.g.ub;
            p = obj.p.sym; p_val = obj.p.val;

            G = obj.G.sym; lbG = zeros(size(G)); ubG = inf*ones(size(G));
            H = obj.H.sym; lbH = zeros(size(H)); ubH = inf*ones(size(H));

            lift_G = SX.sym('lift_G', size(G));
            lift_H = SX.sym('lift_H', size(H));

            ubG(I1) = 0;
            ubH(I2) = 0;
            g = [g;lift_G - G;lift_H - H]; lbg = [lbg;zeros(size(G));zeros(size(H))]; ubg = [ubg;zeros(size(G));zeros(size(H))];
            w = [w;lift_G;lift_H]; lbw = [lbw;lbG;lbH]; ubw = [ubw;ubG;ubH];
            
            bnlp.x = w;
            bnlp.g = g;
            bnlp.p = p;
            bnlp.f = obj.f;

            solver_initialization.lbx = lbw;
            solver_initialization.ubx = ubw;
            solver_initialization.x0 = [obj.w.res;obj.G.eval;obj.H.eval];
            solver_initialization.lbg = lbg;
            solver_initialization.ubg = ubg;
            solver_initialization.p = p_val;
        end

        function [I1, I2] = extract_active_set_from_result(obj, tol)
            arguments
                obj
                tol=1e-4;
            end
            I1 = find(obj.G.eval <= tol);
            I2 = find(obj.G.eval > tol);
        end

        function [nlp,metadata] = get_corresponding_nlp(obj, psi_fun, opts) % TODO(@anton) steering strategy here??
        % Method which generates a relaxed problem with one of the three steering strategies:
        % ($\ell_1$-penalty, $\ell_\infty$-penalty and direct) 
        % and the passed function $\psi(H,G,\sigma)$.
        % TODO (@anton) 
            arguments
                obj
                psi_fun casadi.Function
                opts % TODO(@anton) Maybe this should be a reduced options set that lives in VDX
            end
            import casadi.*
            metadata = struct;
            metadata.ind_map_G = [];
            metadata.ind_map_H = [];
            obj.finalize_assignments();
            casadi_symbolic_mode = obj.w.casadi_type;
            nlp = vdx.Problem('casadi_type', casadi_symbolic_mode);

            % Copy over w, g, p, and f.
            nlp.w = copy(obj.w);
            nlp.w.problem = nlp;
            nlp.g = copy(obj.g);
            nlp.g.problem = nlp;
            nlp.p = copy(obj.p);
            nlp.p.problem = nlp;
            nlp.f = obj.f;

            g_vars = nlp.g.get_vars();
            pre_sort_g_ind = cellfun(@(x) [x.indices{:}], g_vars, "UniformOutput", false);
            pre_sort_g_ind = [pre_sort_g_ind{:}];
            w_vars = nlp.w.get_vars();
            pre_sort_w_ind = cellfun(@(x) [x.indices{:}], w_vars, "UniformOutput", false);
            pre_sort_w_ind = [pre_sort_w_ind{:}];
            p_vars = nlp.p.get_vars();
            pre_sort_p_ind = cellfun(@(x) [x.indices{:}], p_vars, "UniformOutput", false);
            pre_sort_p_ind = [pre_sort_p_ind{:}];

            nlp.p.sigma_p = {{'sigma_p', 1}, opts.sigma_0};

            switch opts.homotopy_steering_strategy
              case HomotopySteeringStrategy.DIRECT
                % Nothing to be done here.
              case HomotopySteeringStrategy.ELL_INF
                % adding a scalar elastic variable to nlp.w which augments the original mpcc.w
                nlp.w.s_elastic = {{'s_elastic', 1}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                if opts.decreasing_s_elastic_upper_bound
                    nlp.g.s_ub = {nlp.w.s_elastic() - nlp.p.sigma_p(), -inf, 0};
                end
                if opts.objective_scaling_direct
                    nlp.f = nlp.f + (1/nlp.p.sigma_p())*nlp.w.s_elastic(); % penalize the elastic more and more with decreasing sigma_p
                else
                    nlp.f = nlp.p.sigma_p()*nlp.f + nlp.w.s_elastic(); % reduce the weight of the initial objective with decreasing sigma_p
                end
              case HomotopySteeringStrategy.ELL_1
                % Ell 1 steering strategy slacks are added during creation of constraints
            end
            % Get all the vdx.Variable for the complementarities
            comp_var_names = obj.G.get_var_names();
            sum_elastic = 0;
            for ii=1:numel(comp_var_names)
                name = comp_var_names{ii};
                depth = obj.G.(name).depth; % how many indices does this variable need?
                s_elastic_name = ['s_elastic_' char(name)];

                % If a 0 dimensional variable handle it specially.
                if depth == 0
                    % Get the corresponding G and H expressions
                    G_curr = obj.G.(name)();
                    H_curr = obj.H.(name)();
                    % select relaxation slacks/parameters
                    switch opts.homotopy_steering_strategy
                      case HomotopySteeringStrategy.DIRECT
                        % nlp.p.sigma_p(): sigma is a parameter/variable that has no indices
                        sigma = nlp.p.sigma_p(); 
                      case HomotopySteeringStrategy.ELL_INF
                        sigma = nlp.w.s_elastic(); % Here s_elastic takes the role of sigma in direct, and sigma_p is used to define a penalty parameter for the elastic variable s_elastic
                      case HomotopySteeringStrategy.ELL_1
                        % adding elastic variables to nlp.w which augments the original mpcc.w
                        % Remark: ELL_1 with s_elastic is equivalent to usually Ell_1 penality approach, but this indirect way helps
                        % to add some constraint on s_elastic (which  avoids unbounded problems somtimes, and it can also improve convergence)
                        nlp.w.(s_elastic_name) = {{s_elastic_name, size(G_curr, 1)}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                        if opts.decreasing_s_elastic_upper_bound
                            % TODO(@anton) this creates a performance bottleneck
                            nlp.g.([s_elastic_name '_ub']) = {nlp.w.(s_elastic_name)() - nlp.p.sigma_p(), -inf, 0};
                        end

                        % TODO(@anton) this creates a performance bottleneck
                        sigma = nlp.w.(s_elastic_name)();
                        sum_elastic = sum_elastic + sum1(sigma);
                    end

                    if opts.lift_complementarities || ~opts.assume_lower_bounds
                        % TODO(@anton) this is in an if for performance reasons not to call find_nonscalar many times
                        %              however it is likely that this can be done in 1 shot?
                        [ind_scalar_G,ind_nonscalar_G, ind_map_G] = find_nonscalar(G_curr, nlp.w.sym);
                        [ind_scalar_H,ind_nonscalar_H, ind_map_H] = find_nonscalar(H_curr, nlp.w.sym);
                        
                        metadata.ind_map_G = [metadata.ind_map_G, ind_map_G];
                        metadata.ind_map_H = [metadata.ind_map_H, ind_map_H];
                    end
                    
                    if opts.lift_complementarities
                        nlp.w.([name '_G_lift']) = {{'G', length(ind_nonscalar_G)}, 0, inf};
                        G_lift = G_curr(ind_nonscalar_G);
                        G_curr(ind_nonscalar_G) = nlp.w.([name '_G_lift'])();

                        nlp.w.([name '_H_lift']) = {{'H', length(ind_nonscalar_H)}, 0, inf};
                        H_lift = H_curr(ind_nonscalar_H);
                        H_curr(ind_nonscalar_H) = nlp.w.([name '_H_lift'])();
                        
                        nlp.g.([name '_G_lift']) = {nlp.w.([name '_G_lift'])()-G_lift};
                        nlp.g.([name '_H_lift']) = {nlp.w.([name '_H_lift'])()-H_lift};    
                    end

                    % If the variables are not already lower bounded we add lower bounds to the nonscalar entries
                    % in the generic constraints and as box constraints to the scalar entries.
                    if ~opts.assume_lower_bounds % Lower bounds on G, H, not already present in MPCC
                        if ~opts.lift_complementarities
                            if ~isempty(ind_nonscalar_G)
                                nlp.g.([name '_G_lb']) = {G_curr(ind_nonscalar_G), 0, inf};
                            end
                            if ~isempty(ind_nonscalar_H)
                                nlp.g.([name '_H_lb']) = {H_curr(ind_nonscalar_H), 0, inf};
                            end
                        end
                        lbw = nlp.w.lb;
                        lbw(ind_map_H) = 0;
                        lbw(ind_map_G) = 0;
                        nlp.w.lb = lbw;
                    end

                    % Build expression and add correct bounds.
                    g_comp_expr = psi_fun(G_curr, H_curr, sigma);
                    [lb, ub, g_comp_expr] = generate_mpcc_relaxation_bounds(g_comp_expr, obj.opts.relaxation_strategy);
                    nlp.g.([name '_relaxed']) = {g_comp_expr,lb,ub};
                else % Do the same behavior as before excep this time for each var(i,j,k,...) index for each variable 'var'.
                     % Get indices that we will need to get all the casadi vars for the vdx.Variable
                    indices = {};
                    for len=size(obj.G.(name).indices)
                        indices = horzcat(indices, {1:len});
                    end
                    indices = indices(1:depth);
                    % We subtract 1 to get the 0 indexing correct :)
                    inorderlst = all_combinations(indices{:})-1;

                    % Transfer each element of the vdx.Variable individually.
                    for ii=1:size(inorderlst,1)
                        curr = num2cell(inorderlst(ii,:));
                        G_curr = obj.G.(name)(curr{:});
                        H_curr = obj.H.(name)(curr{:});
                        if any(size(G_curr) == 0)
                            continue
                        end
                        % select relaxation slacks/parameters
                        switch opts.homotopy_steering_strategy
                          case HomotopySteeringStrategy.DIRECT
                            % nlp.p.sigma_p(): sigma is a parameter/variable that has no indices
                            sigma = nlp.p.sigma_p(); 
                          case HomotopySteeringStrategy.ELL_INF
                            sigma = nlp.w.s_elastic(); % Here s_elastic takes the role of sigma in direct, and sigma_p is used to define a penalty parameter for the elastic variable s_elastic
                          case HomotopySteeringStrategy.ELL_1
                            % adding elastic variables to nlp.w which augments the original mpcc.w
                            % Remark: ELL_1 with s_elastic is equivalent to usually Ell_1 penality approach, but this indirect way helps
                            % to add some constraint on s_elastic (which  avoids unbounded problems somtimes, and it can also improve convergence)
                            nlp.w.(s_elastic_name)(curr{:}) = {{s_elastic_name, size(G_curr, 1)}, opts.s_elastic_min, opts.s_elastic_max, opts.s_elastic_0};
                            if opts.decreasing_s_elastic_upper_bound
                                % TODO(@anton) this creates a performance bottleneck
                                nlp.g.([s_elastic_name '_ub'])(curr{:}) = {nlp.w.(s_elastic_name)(curr{:}) - nlp.p.sigma_p(), -inf, 0};
                            end

                            % TODO(@anton) this creates a performance bottleneck
                            sigma = nlp.w.(s_elastic_name)(curr{:});
                            sum_elastic = sum_elastic + sum1(sigma);
                        end

                        if opts.lift_complementarities || ~opts.assume_lower_bounds
                            % TODO(@anton) this is in an if for performance reasons not to call find_nonscalar many times
                            %              however it is likely that this can be done in 1 shot?
                            [ind_scalar_G,ind_nonscalar_G, ind_map_G] = find_nonscalar(G_curr, nlp.w.sym);
                            [ind_scalar_H,ind_nonscalar_H, ind_map_H] = find_nonscalar(H_curr, nlp.w.sym);

                            metadata.ind_map_G = [metadata.ind_map_G, ind_map_G];
                            metadata.ind_map_H = [metadata.ind_map_H, ind_map_H];
                        end
                        
                        if opts.lift_complementarities
                            nlp.w.([name '_G_lift'])(curr{:}) = {{'G', length(ind_nonscalar_G)}, 0, inf};
                            G_lift = G_curr(ind_nonscalar_G);
                            G_curr(ind_nonscalar_G) = nlp.w.([name '_G_lift'])(curr{:});

                            nlp.w.([name '_H_lift'])(curr{:}) = {{'H', length(ind_nonscalar_H)}, 0, inf};
                            H_lift = H_curr(ind_nonscalar_H);
                            H_curr(ind_nonscalar_H) = nlp.w.([name '_H_lift'])(curr{:});
                            
                            nlp.g.([name '_G_lift'])(curr{:}) = {nlp.w.([name '_G_lift'])(curr{:})-G_lift};
                            nlp.g.([name '_H_lift'])(curr{:}) = {nlp.w.([name '_H_lift'])(curr{:})-H_lift};    
                        end

                        if ~opts.assume_lower_bounds % Lower bounds on G, H, not already present in MPCC
                            if ~opts.lift_complementarities
                                if ~isempty(ind_nonscalar_G)
                                    nlp.g.([name '_G_lb'])(curr{:}) = {G_curr(ind_nonscalar_G), 0, inf};
                                end
                                if ~isempty(ind_nonscalar_H)
                                    nlp.g.([name '_H_lb'])(curr{:}) = {H_curr(ind_nonscalar_H), 0, inf};
                                end
                            end
                            lbw = nlp.w.lb;
                            lbw(ind_map_H) = 0;
                            lbw(ind_map_G) = 0;
                            nlp.w.lb = lbw;
                        end
                        
                        g_comp_expr = psi_fun(G_curr, H_curr, sigma);
                        [lb, ub, g_comp_expr] = generate_mpcc_relaxation_bounds(g_comp_expr, opts.relaxation_strategy);
                        nlp.g.(name)(curr{:}) = {g_comp_expr,lb,ub};
                    end
                end
            end

            % add ell_1 cost if necessary
            if opts.objective_scaling_direct
                nlp.f = nlp.f + (1/nlp.p.sigma_p())*sum_elastic;
            else
                nlp.f = nlp.p.sigma_p()*nlp.f + sum_elastic;
            end

            % Build maps to recover original order from the nlp.
            nlp.g.sort_by_index;
            w_map = nlp.w.sort_by_index;
            if ~isempty(metadata.ind_map_G)
                [metadata.ind_map_G,~] = find(w_map == metadata.ind_map_G);
            end
            if ~isempty(metadata.ind_map_H)
                [metadata.ind_map_H,~] = find(w_map == metadata.ind_map_H);
            end
            post_sort_g_ind = cellfun(@(x) [x.indices{:}], g_vars, "UniformOutput", false);
            post_sort_g_ind = [post_sort_g_ind{:}];
            metadata.ind_map_g.mpcc = pre_sort_g_ind;
            metadata.ind_map_g.nlp = post_sort_g_ind;
            post_sort_w_ind = cellfun(@(x) [x.indices{:}], w_vars, "UniformOutput", false);
            post_sort_w_ind = [post_sort_w_ind{:}];
            metadata.ind_map_w.mpcc = pre_sort_w_ind;
            metadata.ind_map_w.nlp = post_sort_w_ind;
            post_sort_p_ind = cellfun(@(x) [x.indices{:}], p_vars, "UniformOutput", false);
            post_sort_p_ind = [post_sort_p_ind{:}];
            metadata.ind_map_p.mpcc = pre_sort_p_ind;
            metadata.ind_map_p.nlp = post_sort_p_ind;

            % Build required casadi functions for recovering mpcc data from the nlp.
            G_fun = Function('G', {nlp.w.sym, nlp.p.sym}, {obj.G.sym});
            H_fun = Function('H', {nlp.w.sym, nlp.p.sym}, {obj.H.sym});
            metadata.G_fun = G_fun;
            metadata.H_fun = H_fun;
            comp_res_fun = Function('comp_res', {nlp.w.sym, nlp.p.sym}, {mmax(obj.G.sym.*obj.H.sym)});
            metadata.comp_res_fun = comp_res_fun;
            f_mpcc_fun = Function('f_mpcc', {nlp.w.sym, nlp.p.sym}, {obj.f});
            metadata.f_mpcc_fun = f_mpcc_fun;
            w_mpcc_fun = Function('w_mpcc', {nlp.w.sym}, {obj.w.sym});
            metadata.w_mpcc_fun = w_mpcc_fun;
            g_mpcc_fun = Function('g_mpcc', {nlp.w.sym, nlp.p.sym}, {obj.g.sym});
            metadata.g_mpcc_fun = g_mpcc_fun;
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

function X = all_combinations(varargin)
    numSets = length(varargin);
    for i=1:numSets,
        thisSet = sort(varargin{i});
        if ~isequal(prod(size(thisSet)),length(thisSet)),
            nosnoc.error('combinations_not_vectors', 'All inputs must be vectors.')
        end
        if ~isnumeric(thisSet),
            nosnoc.error('cominations_not_numeric','All inputs must be numeric.')
        end
        sizeThisSet(i) = length(thisSet);
        varargin{i} = thisSet;
    end
    X = zeros(prod(sizeThisSet),numSets);
    for i=1:size(X,1)
        ixVect = cell(length(sizeThisSet),1);
        [ixVect{:}] = ind2sub(flip(sizeThisSet),i);
        ixVect = flip([ixVect{:}]);
        vect = zeros(1, numSets);
        for jj=1:numSets
            vect(jj) = varargin{jj}(ixVect(jj));
        end
        X(i,:) = vect;
    end
end
