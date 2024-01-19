% Solve an example QP via FESD applied to the vector flow.
% TODO(@anton) I want to see how horrible this is on large QPs.
Q = -2*eye(2);
c = [2;2];
A = [eye(2);-1,-1];
b = [0;0;-4];
x0 = [1.1;1.05];
results = qp_solver(Q,c,A,b,x0)
