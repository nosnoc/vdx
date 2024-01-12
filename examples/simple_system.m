clear all
close all
import casadi.*
import vdx.*

T = 2.0;
N_sim = 100;
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
data.c = [x(2);-x(2) - (x(1)+0.5)^2 + 1];
data.f_x = [x(2); -x(1)];
data.f_q = 0;
data.f_q_T = 0;

data.T = t_step;
data.N_stages = 1;
data.N_fe = 3;
data.n_s = 2;
data.irk_scheme = 'radau';

prob = InclusionProblem(data, struct);

prob.generate_constraints();

prob.create_solver(struct);

x_res = data.x0;
h_res = [];
x_curr = data.x0;
for step=1:N_sim
    prob.w.x(0,0,0).init = x_curr;
    prob.w.x(0,0,0).lb = x_curr;
    prob.w.x(0,0,0).ub = x_curr;
    homotopy(prob)
    x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
    x_sim = prob.w.x(1,:,data.n_s).res;
    x_sim = [x_sim{:}];
    h_sim = prob.w.h(1,:).res;
    h_sim = [h_sim{:}];
    x_res = [x_res,x_sim];
    h_res = [h_res,h_sim];
end

t_res = [0,cumsum(h_res)];

figure
hold on
plot(t_res, x_res(1,:))
plot(t_res, x_res(2,:))
hold off
figure
stairs(t_res(1:end-1),h_res)
figure
hold on
x1 = -1.7:0.001:0.7;
x2 = -(x1+0.5).^2 + 1;
plot(x_res(1,:), x_res(2,:))
plot(x1,x2,'r--')
yline(0, 'r--')
hold off
