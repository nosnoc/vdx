classdef InclusionProblem < vdx.Problem
    properties (Access=public)
        data
        opts
        eta_vec
        eta_fun
    end
    
    methods (Access=public)
        function obj = InclusionProblem(data, opts)
            obj@vdx.Problem(); % call parent constructor to prepare
            obj.data = data;
            obj.opts = opts;

            if ~isfield(obj.opts, 'use_fesd')
                obj.opts.use_fesd = true;
            end
            if ~isfield(obj.opts,'comp_scale')
                obj.opts.comp_scale = 1;
            end
            % get dimensions
            n_x = length(data.x);
            n_u = length(data.u);
            n_c = length(data.c);

            % if opts.elastic_ell_inf
            %     obj.w.s_elastic(1) = {{'s_elastic',1},0,inf,0};
            % end
            obj.p.sigma(1) = {{'sigma',1},0,inf,0};
            obj.p.gamma_h(1) = {{'gamma_h',1},0,inf,0.1};
            obj.p.T(1) = {{'T',1},0,inf,data.T};


            % other derived values
            t_stage = data.T/data.N_stages;
            h0 = t_stage/data.N_fe;
            
            obj.w.x(0,0,data.n_s) = {{['x_0'], n_x}};
            for ii=1:data.N_stages
                obj.w.u(ii) = {{['u_' num2str(ii)], n_u}, data.lbu, data.ubu, data.u0};
                for jj=1:data.N_fe
                    if obj.opts.use_fesd
                        obj.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, 0, 2*h0, h0};
                    end
                    for kk=1:data.n_s
                        obj.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_x}, data.lbx, data.ubx, data.x0};
                        obj.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c},0,inf};
                    end
                end
            end
        end

        function generate_constraints(obj)
            import casadi.*
            [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.data.n_s, obj.data.irk_scheme);
            n_x = length(obj.data.x);
            n_u = length(obj.data.u);
            n_c = length(obj.data.c);

            % other derived values
            if obj.opts.use_fesd
                t_stage = obj.p.T(1)/obj.data.N_stages;
                h0 = obj.p.T(1).init/(obj.data.N_stages*obj.data.N_fe);
            else
                h0 = obj.p.T(1).init/(obj.data.N_stages*obj.data.N_fe);
            end
            
            % Define functions from obj.data
            lambda = SX.sym('lambda', n_c);

            nabla_c = obj.data.c.jacobian(obj.data.x)';

            f_x = obj.data.f_x + nabla_c*lambda;

            f_x_fun = Function('f_x_fun', {obj.data.x,obj.data.u,lambda}, {f_x});
            f_q_fun = Function('q_fun', {obj.data.x,obj.data.u}, {obj.data.f_q});
            f_q_T_fun = Function('q_fun', {obj.data.x}, {obj.data.f_q_T});
            c_fun = Function('c_fun', {obj.data.x}, {obj.data.c});
            x_prev = obj.w.x(0,0,obj.data.n_s);
            for ii=1:obj.data.N_stages
                ui = obj.w.u(ii);
                sum_h = 0;
                for jj=1:obj.data.N_fe
                    if obj.opts.use_fesd
                        h = obj.w.h(ii,jj);
                        sum_h = sum_h + h;
                    else
                        h = h0;
                    end
                    for kk=1:obj.data.n_s
                        x_ijk = obj.w.x(ii,jj,kk);
                        lambda_ijk = obj.w.lambda(ii,jj,kk);
                        fj = f_x_fun(x_ijk,ui,lambda_ijk);
                        qj = f_q_fun(x_ijk,ui);
                        xk = C(1, kk+1) * x_prev;
                        for rr=1:obj.data.n_s
                            x_ijr = obj.w.x(ii,jj,rr);
                            xk = xk + C(rr+1, kk+1) * x_ijr;
                        end
                        obj.g.dynamics(ii,jj,kk) = {h * fj - xk};
                        % also add non-negativity constraint on c
                        obj.g.c_nonnegative(ii,jj,kk) = {c_fun(x_ijk), 0, inf};
                        % also integrate the objective
                        obj.f = obj.f + B(kk+1)*h*qj;
                    end
                    x_prev = obj.w.x(ii,jj,obj.data.n_s);
                end
                if obj.opts.use_fesd
                    obj.g.sum_h(ii) = {t_stage-sum_h};
                end
            end

            % Terminal cost
            obj.f = obj.f + f_q_T_fun(obj.w.x(ii,jj,kk));

            % Terminal constraint
            if isfield(obj.data, 'g_T')
                g_T_fun = Function('g_T_fun', {obj.data.x}, {obj.data.g_T});
                obj.g.terminal(0) = {g_T_fun(obj.w.x(ii,jj,kk))}; % TODO(@anton) assume equality for now
            end

            % C indicator implementation.
            if 0 %TODO(@anton) this seems to work. remove once sure
                x = obj.w.x(0,0,obj.data.n_s);
                cx = c_fun(x);
                obj.w.c_ind(0,0,obj.data.n_s) = {{'c_ind_0', n_c}, 0, 1, 0};
                c_ind = obj.w.c_ind(0,0,obj.data.n_s);
                obj.g.c_ind_comp(0,0,obj.data.n_s) = {cx.*c_ind - obj.p.sigma(1), -inf, 0};
                obj.f = obj.f + sum1(1e-2*obj.p.gamma_h(1)*(c_ind-1).^2);
                for ii=1:obj.data.N_stages
                    for jj=1:obj.data.N_fe
                        for kk=1:obj.data.n_s
                            x = obj.w.x(ii,jj,kk);
                            cx = c_fun(x);
                            obj.w.c_ind(ii,jj,kk) = {{['c_ind_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c}, 0, inf, 0};
                            c_ind = obj.w.c_ind(ii,jj,kk);
                            obj.g.c_ind_comp(ii,jj,kk) = {cx.*c_ind - obj.p.sigma(1), -inf, 0};
                            obj.f = obj.f + sum1(1e-2*obj.p.gamma_h(1)*(c_ind-1).^2);
                        end
                    end
                end
            end

            
            % Do Cross-Complementarity
            x_prev = obj.w.x(0,0,obj.data.n_s);
            G = [];
            H = [];
            for ii=1:obj.data.N_stages
                for jj=1:obj.data.N_fe
                    if obj.opts.use_fesd
                        Gij = c_fun(x_prev);
                        Hij = 0;
                        for kk=1:obj.data.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            Gij = Gij + c_fun(x_ijk);
                            Hij = Hij + lambda_ijk;
                        end
                        G = [G;Gij];
                        H = [H;Hij];
                        obj.g.complementarity(ii,jj) = {obj.opts.comp_scale*(Gij.*Hij - obj.p.sigma(1)), -inf, 0};
                        x_prev = obj.w.x(ii,jj,obj.data.n_s);
                    else
                        Gij = [];
                        Hij = [];
                        for kk=1:obj.data.n_s
                            x_ijk = obj.w.x(ii,jj,kk);
                            lambda_ijk = obj.w.lambda(ii,jj,kk);
                            Gij = [Gij;c_fun(x_ijk)];
                            Hij = [Hij;lambda_ijk];
                        end
                        G = [G;Gij];
                        H = [H;Hij];
                        obj.g.complementarity(ii,jj) = {obj.opts.comp_scale*(Gij.*Hij - obj.p.sigma(1)), -inf, 0};
                    end
                end
            end

            
            % Do Step equilibration (only heuristic for now)
            if ~isfield(obj.opts, 'step_eq')
                obj.opts.step_eq = 'heuristic_mean';
            end
            if ~obj.opts.use_fesd
                obj.opts.step_eq = 'none';
            end
            switch obj.opts.step_eq
              case 'heuristic_mean'
                for ii=1:obj.data.N_stages
                    for jj=1:obj.data.N_fe
                        obj.f = obj.f + obj.p.gamma_h(1)*(h0-obj.w.h(ii,jj))^2;
                    end
                end
              case 'direct'
                eta_vec = [];
                for ii=1:obj.data.N_stages
                    if ii > 1
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii-1,jj,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii-1,jj,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,1,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,1,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,1) - obj.w.h(ii-1,jj);
                        obj.g.step_equilibration(ii,1) = {eta*delta_h, -1e-6, 1e-6};
                    end
                    for jj=2:obj.data.N_fe
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        obj.g.step_equilibration(ii,jj) = {eta*delta_h, -1e-6, 1e-6};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});
              case 'direct_homotopy'
                eta_vec = [];
                for ii=1:obj.data.N_stages
                    if ii > 1
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii-1,jj,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii-1,jj,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,1,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,1,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,1) - obj.w.h(ii-1,jj);
                        homotopy_eq = [eta*delta_h - obj.p.sigma(1);eta*delta_h + obj.p.sigma(1)];
                        obj.g.step_equilibration(ii,1) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                    for jj=2:obj.data.N_fe
                        sigma_c_B = 0;
                        sigma_lam_B = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_B = sigma_c_B + c_fun(obj.w.x(ii,jj-1,kk));
                            sigma_lam_B = sigma_lam_B + obj.w.lambda(ii,jj-1,kk);
                        end
                        sigma_c_F = 0;
                        sigma_lam_F = 0;
                        for kk=1:obj.data.n_s
                            sigma_c_F = sigma_c_F + c_fun(obj.w.x(ii,jj,kk));
                            sigma_lam_F = sigma_lam_F + obj.w.lambda(ii,jj,kk);
                        end

                        pi_c = sigma_c_B .* sigma_c_F;
                        pi_lam = sigma_lam_B .* sigma_lam_F;
                        nu = pi_c + pi_lam;
                        eta = 1;
                        for jjj=1:length(nu)
                            eta = eta*nu(jjj);
                        end
                        eta_vec = [eta_vec;eta];
                        obj.eta_vec = eta_vec;
                        delta_h = obj.w.h(ii,jj) - obj.w.h(ii,jj-1);
                        homotopy_eq = [eta*delta_h - obj.p.sigma(1);eta*delta_h + obj.p.sigma(1)];
                        obj.g.step_equilibration(ii,jj) = {homotopy_eq, [-inf;0], [0;inf]};
                    end
                end
                obj.eta_fun = Function('eta_fun', {obj.w.w}, {eta_vec});

              case 'direct_fix_pathological'
                %TODO(@anton)
              case 'none'
                % Do nothing
            end
        end
    end
end
