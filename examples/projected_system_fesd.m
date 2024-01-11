clear all
close all
import casadi.*
import vdx.*

%% Define (uncontrolled for now) projected system
x = SX.sym('x', 2);
c = [x(2);-x(2) - x(1)^2 + 1];
f = [x(2); -x(1)];

%% Define the needed functions
n_x = length(x);
n_c = length(c);

lambda = SX.sym('lambda', n_c);

nabla_c = c.jacobian(x);

f_x = f + nabla_c*lambda;

f_x_fun = Function('f_x_fun', {x,lambda}, {f_x});
c_fun = Function('c_fun', {x}, {c});

%% FESD problem data
N_stages = 1;
N_fe = 3;
n_s = 2;

T = 0.1;
h0 = T/N_fe;
irk_scheme = 'radau';
[B, C, D, tau_root] = generate_butcher_tableu_integral(n_s, irk_scheme);

prob = Problem();

%% Create variables for problem
prob.w.x(0,0,0) = {{['x_0'], n_x}};
for ii=1:N_stages
    for jj=1:N_fe
        prob.w.h(ii,jj) = {{['h_' num2str(ii) '_' num2str(jj)], 1}, 0, 0.1};
        for kk=1:n_s
            prob.w.x(ii,jj,kk) = {{['x_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_x}};
            prob.w.lambda(ii,jj,kk) = {{['lambda_' num2str(ii) '_' num2str(jj) '_' num2str(kk)], n_c},0,inf};
        end
    end
end

prob.p.sigma(1) = {SX.sym('sigma'),0,inf,0};

%% Dyanamics
x_prev = prob.w.x(0,0,0);
for ii=1:N_stages
    for jj=1:N_fe
        for kk=1:n_s
            x_ijk = prob.w.x(ii,jj,kk);
            lambda_ijk = prob.w.lambda(ii,jj,kk);
            fj = f_x_fun(x_ijk,lambda_ijk);
            xk = C(1, kk+1) * x_prev;
            for rr=1:n_s
                x_ijr = prob.w.x(ii,jj,rr);
                xj = xk + C(rr+1, kk+1) * x_ijr;
            end
            h = prob.w.h(ii,jj);
            prob.g.dynamics(ii,jj,kk) = {h * fj - xj};
            % also add non-negativity constraint on c
            prob.g.c_nonnegative(ii,jj,kk) = {c_fun(x_ijk), 0, inf};
        end
        x_prev = prob.w.x(ii,jj,n_s);
    end
end

%% Cross Complementarity
% In this case using FE-FE only because its easy to implement :)

x_prev = prob.w.x(0,0,0);
G = [];
H = [];
for ii=1:N_stages
    for jj=1:N_fe
        Gij = c_fun(x_prev);
        Hij = 0;
        for kk=1:n_s
            x_ijk = prob.w.x(ii,jj,kk);
            lambda_ijk = prob.w.lambda(ii,jj,kk);
            Gij = Gij + c_fun(x_ijk);
            Hij = Hij + lambda_ijk;
        end
        G = [G;Gij];
        H = [H;Hij];
        prob.g.complementarity(ii,jj) = {Gij.*Hij - prob.p.sigma(1), 0, inf};
        x_prev = prob.w.x(ii,jj,n_s);
    end
end


