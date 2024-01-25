clear all
close all
import casadi.*
import vdx.*

T = 1;
R = 3.5;
%% Define projected system
x1 = SX.sym('x1', 2);
x2 = SX.sym('x2', 2);
x = [x1;x2];
x_target = [10;0;10;10];
data.x = x;
data.lbx = [-inf;-inf;-inf;-inf];
data.ubx = [inf;inf;inf;inf];
data.x0 = [-10;0;10;0];
u1 = SX.sym('u1', 2);
u2 = SX.sym('u2', 2);
data.u = [u1;u2];
data.lbu = [-100/sqrt(2);-100/sqrt(2);0;0];
%data.lbu = [-100/sqrt(2);-100/sqrt(2);-60/sqrt(2);-60/sqrt(2)];
%data.ubu = [100/sqrt(2);100/sqrt(2);60/sqrt(2);60/sqrt(2)];
data.ubu = [100/sqrt(2);100/sqrt(2);0;0];
data.u0 = [0;0;0;0];
data.c = [norm_2(x2-x1)-2*R];
data.f_x = [u1;0;0];

% costs
data.f_q = 0.001*norm_2(data.u)^2;
data.f_q_T = 10000*0.5*(norm_2(x2-x_target(3:4))^2);%0.5*(norm_2(x)^2);

data.T = T;
data.N_stages = 25;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

opts.step_eq = 'direct';

prob = InclusionProblem(data, opts);

prob.generate_constraints();

%% create solver
default_tol = 1e-12;

%opts_casadi_nlp.ipopt.print_level = 1;
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

%% Do homotopy
prob.w.x(0,0,data.n_s).init = data.x0;
prob.w.x(0,0,data.n_s).lb = data.x0;
prob.w.x(0,0,data.n_s).ub = data.x0;
homotopy(prob);
%% plot
x_res = prob.w.x(0:data.N_stages,0:data.N_fe,data.n_s).res';
x_res = [x_res{:}];
u_res = prob.w.u(1:data.N_stages).res';
u_res = [u_res{:}];
h_res = prob.w.h(:,:).res';
h_res = [h_res{:}];
t_res = [0,cumsum(h_res)]
plot_discs(h_res,x_res,[3.5,3.5], ["circle", "circle"])
