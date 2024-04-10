import casadi.*
T = 10; % Time horizon
N = 20; % number of control intervals

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
xnext = SX.sym('y',2);
u = SX.sym('u');
p = SX.sym('p')

% Model equations
xdot = [(1-x2^2)*x1 - x2 + u; p*x1];

% Objective term
L = x1^2 + x2^2 + u^2;


% Fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = T/N/M;
f = Function('f', {x, u, p}, {xdot, L});
X0 = MX.sym('X0', 2);
X1 = MX.sym('X1', 2);
U = MX.sym('U');
P = MX.sym('P');
X = X0;
Q = 0;
for j=1:M
    [k1, k1_q] = f(X, U, P);
    [k2, k2_q] = f(X + DT/2 * k1, U, P);
    [k3, k3_q] = f(X + DT/2 * k2, U, P);
    [k4, k4_q] = f(X + DT * k3, U, P);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = Function('F', {X0, U, P}, {X, Q}, {'x0','u', 'p'}, {'xf', 'qf'});
dyn = Function('dyn', {X1, X0, U, P}, {X-X1});
qf = Function('qf', {X0, U, P}, {Q});
g_path = Function('g_path', {x}, {x1+x2+2});

warning off vdx:indexing:dot_reference_returns_vdx_var % Turn off warnings for advanced features of vdx.
% create a problem
prob = vdx.Problem('casadi_type', 'MX');

prob.p.p = {P, -inf, inf, 1};
prob.w.x(0) = {{'X_0', 2},[0;1], [0;1], [0;1]};
prob.w.x(1:N) = {{'X', 2}, [-0.25;-inf], [inf;inf], [0;0]};
prob.w.u(1:N) = {{'U', 1},-1,1,0};
prob.w.add_variable_group('stage_vars',{prob.w.x, prob.w.u});
prob.w.add_variable_group('dynamics_args',{prob.w.x, prob.w.x, prob.w.u}, {[],@vdx.indexing.previous,[]});
prob.g.path(0:N) = {{g_path, {prob.w.x}},0,inf};
prob.g.dynamics(1:N) = {{dyn, {prob.w.x, prob.w.x, prob.w.u, prob.p.p}, {[],@vdx.indexing.previous,[],[]}}};
Xk = prob.w.x(0);
% Formulate the NLP
for k=1:N
    % New NLP variable for the control
    Uk = prob.w.u(k);

    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'u', Uk, 'p', prob.p.p(0));
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

warning on vdx:indexing:dot_reference_returns_vdx_var % Turn back on warnings for advanced features of vdx.
