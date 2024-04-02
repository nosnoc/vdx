import casadi.*
T = 10; % Time horizon
N = 20; % number of control intervals

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
u = SX.sym('u');

% Model equations
xdot = [(1-x2^2)*x1 - x2 + u; x1];

% Objective term
L = x1^2 + x2^2 + u^2;

% Continuous time dynamics
f = Function('f', {x, u}, {xdot, L});


% Fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = T/N/M;
f = Function('f', {x, u}, {xdot, L});
X0 = MX.sym('X0', 2);
X1 = MX.sym('X1', 2);
U = MX.sym('U');
X = X0;
Q = 0;
for j=1:M
    [k1, k1_q] = f(X, U);
    [k2, k2_q] = f(X + DT/2 * k1, U);
    [k3, k3_q] = f(X + DT/2 * k2, U);
    [k4, k4_q] = f(X + DT * k3, U);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});
dyn = Function('dyn', {X0, X1, U}, {X-X1});
qf = Function('qf', {X0, U}, {Q});
g_path = Function('g_path', {x}, {x1+x2+2});

% create a problem
prob = vdx.Problem('casadi_type', 'MX');

% "Lift" initial conditions
prob.w.x(0) = {{'X_0', 2},[0;1], [0;1], [0;1]};
prob.w.x(1:N) = {{'X', 2}, [-0.25;-inf], [inf;inf], [0;0]};
prob.w.u(1:N) = {{'U', 1},-1,1,0};
prob.w.add_variable_group('stage_vars',["x", "u"]);
prob.w.add_variable_group('dynamics_args',["x", "x", "u"], {[],@vdx.indexing.previous,[]});
prob.g.path(0:N) = {{g_path, ["x"], prob.w, {}},0,inf};
prob.g.dynamics(1:N) = {{dyn, ["x", "x", "u"], prob.w, {[],@vdx.indexing.previous,[]}},0,inf};
Xk = prob.w.x(0);
% Formulate the NLP
for k=1:N
    % New NLP variable for the control
    Uk = prob.w.u(k);

    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'p', Uk);
    prob.f=prob.f+Fk.qf;

    % New NLP variable for state at end of interval
    Xk = prob.w.x(k);
end

% re-sort for that sweet sweet block diagonal structure
prob.w.sort_by_index();
prob.g.sort_by_index();

% Create an NLP solver
casadi_opts = struct;
prob.create_solver(casadi_opts);

% Solve the NLP
prob.solve();

% Plot the solution
x_res = prob.w.x.res;
u_res = prob.w.u.res;
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x_res(1,:), '--')
plot(tgrid, x_res(2,:), '-')
stairs(tgrid, [u_res,nan], '-.')
xlabel('t')
legend('x1','x2','u')
