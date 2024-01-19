% Solve an example QP via FESD applied to the vector flow.
Q = eye(2);
c = [2;2];
A = eye(2);
b = [0;0];
x0 = [0.1;0.1];
results = qp_solver(Q,c,A,b,x0)
