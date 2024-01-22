classdef InclusionProblem < vdx.Problem
    properties (Access=public)
        data
        opts
    end
    
    methods (Access=public)
        function obj = InclusionProblem(data, opts)
            obj@vdx.Problem(); % call parent constructor to prepare
            obj.data = data;
            obj.opts = opts;
            
            % get dimensions
            n_x = length(data.x);
            n_u = length(data.u);
            n_c = length(data.c);

            % other derived values
            t_stage = data.T/data.N_stages;
            h0 = t_stage/data.N_fe;
            
            obj.w.x(0,0,data.n_s) = {{['x_0'], n_x}};
            for ii=1:data.N_stages
                obj.w.u(ii) = {{['u_' num2str(ii)], n_u}, data.lbu, data.ubu, data.u0};
                for jj=1:data.N_fe
                    obj.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, 0, 2*h0, h0};
                    for kk=1:data.n_s
                        obj.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_x}, data.lbx, data.ubx, data.x0};
                        obj.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c},0,inf};
                    end
                end
            end
            obj.p.sigma(1) = {{'sigma',1},0,inf,0};
            obj.p.gamma_h(1) = {{'gamma_h',1},0,inf,0.1};
            obj.p.T(1) = {{'T',1},0,inf,data.T};
        end

        function generate_constraints(obj)
            import casadi.*
            [B, C, D, tau_root] = generate_butcher_tableu_integral(obj.data.n_s, obj.data.irk_scheme);
            n_x = length(obj.data.x);
            n_u = length(obj.data.u);
            n_c = length(obj.data.c);

            % other derived values
            t_stage = obj.p.T(1)/obj.data.N_stages;
            h0 = obj.p.T(1).init/obj.data.N_fe;
            
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
                    h = obj.w.h(ii,jj);
                    sum_h = sum_h + h;
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
                obj.g.sum_h(ii) = {t_stage-sum_h};
            end

            % Terminal cost
            obj.f = obj.f + f_q_T_fun(obj.w.x(ii,jj,kk));

            % Terminal constraint
            if isfield(obj.data, 'g_T')
                g_T_fun = Function('g_T_fun', {obj.data.x}, {obj.data.g_T});
                obj.g.terminal(0) = {g_T_fun(obj.w.x(ii,jj,kk))}; % TODO(@anton) assume equality for now
            end

            % Do Cross-Complementarity
            x_prev = obj.w.x(0,0,obj.data.n_s);
            G = [];
            H = [];
            for ii=1:obj.data.N_stages
                for jj=1:obj.data.N_fe
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
                    obj.g.complementarity(ii,jj) = {Gij.*Hij - obj.p.sigma(1), -inf, 0};
                    x_prev = obj.w.x(ii,jj,obj.data.n_s);
                end
            end

            % Do Step equilibration (only heuristic for now)
            for ii=1:obj.data.N_stages
                for jj=1:obj.data.N_fe
                    obj.f = obj.f + obj.p.gamma_h(1)*(h0-obj.w.h(ii,jj))^2;
                end
            end
        end
    end
end
