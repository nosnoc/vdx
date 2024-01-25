clear all
close all
import casadi.*
import vdx.*

T = 5.0;
N_sim = 20;
t_step = T/N_sim;
%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
data.x = x;
data.lbx = [-inf;-inf];
data.ubx = [inf;inf];
data.x0 = [0;0.5];
data.u = [];
data.lbu = [];
data.ubu = [];
data.u0 = [];
data.c = [x(2)+0.25;-x(2) - (x(1)+0.5)^2 + 1];
data.f_x = [x(2); -x(1)];
data.f_q = 0;
data.f_q_T = 0;

data.T = t_step;
data.N_stages = 1;
data.N_fe = 2;
data.n_s = 1;
data.irk_scheme = 'radau';

opts.step_eq = 'heuristic_mean';
opts.use_fesd = false;

prob = InclusionProblem(data, opts);

prob.generate_constraints();

default_tol = 1e-12;

opts_casadi_nlp.ipopt.print_level = 2;
opts_casadi_nlp.print_time = 0;
opts_casadi_nlp.ipopt.sb = 'yes';
opts_casadi_nlp.verbose = false;
opts_casadi_nlp.ipopt.max_iter = 10000;
opts_casadi_nlp.ipopt.bound_relax_factor = 0;
%opts_casadi_nlp.ipopt.bound_relax_factor = 1e-8;
%opts_casadi_nlp.ipopt.honor_original_bounds = 'yes';
opts_casadi_nlp.ipopt.tol = default_tol;
opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
opts_casadi_nlp.ipopt.dual_inf_tol = default_tol;
opts_casadi_nlp.ipopt.compl_inf_tol = default_tol;
opts_casadi_nlp.ipopt.acceptable_tol = 1e-6;
opts_casadi_nlp.ipopt.mu_strategy = 'adaptive';
opts_casadi_nlp.ipopt.mu_oracle = 'quality-function';
opts_casadi_nlp.ipopt.warm_start_init_point = 'yes';
opts_casadi_nlp.ipopt.warm_start_entire_iterate = 'yes';
opts_casadi_nlp.ipopt.linear_solver = 'ma27';
prob.create_solver(opts_casadi_nlp);

x_res = data.x0;
h_res = [];
x_curr = data.x0;
for step=1:N_sim
    prob.w.x(0,0,data.n_s).init = x_curr;
    prob.w.x(0,0,data.n_s).lb = x_curr;
    prob.w.x(0,0,data.n_s).ub = x_curr;
    success = homotopy(prob);
    if ~success
        disp(['Failure to converge at step=' num2str(step)])
    else
        disp(['step=' num2str(step)])
    end
    x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
    x_sim = prob.w.x(1,:,data.n_s).res;
    x_sim = [x_sim{:}];
    if opts.use_fesd
        h_sim = prob.w.h(1,:).res;
        h_sim = [h_sim{:}];
    else
        h_sim = prob.p.T(1).init/(data.N_fe)*ones(data.N_fe,1)';
    end
    h_res = [h_res,h_sim];
    x_res = [x_res,x_sim];
end

t_res = [0,cumsum(h_res)];

figure
hold on
plot(t_res, x_res(1,:))
plot(t_res, x_res(2,:))
hold off
if opts.use_fesd
    figure
    stairs(t_res(1:end-1),h_res)
end
figure
hold on
x1 = -1.7:0.001:0.7;
x2 = -(x1+0.5).^2 + 1;
plot(x_res(1,:), x_res(2,:))
plot(x1,x2,'r--')
yline(-0.25, 'r--')
hold off
