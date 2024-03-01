%clear all
close all
import vdx.*
N_sim = 10;
ts = 0;
xs = -2.4:0.01:-2;%-1.6:0.1:-1;
x_term = [];

figure
hold on
fplot(@(x) 1/24*(256 + 96*x + 12*x.^2 + x.^3) - (x.^3)/8 + (4+x/2 - 3.5).^2, "DisplayName", 'analytical')
xlim([-2.4 -2])
legend('location', 'north')

for comp_tol=[1,1e-1,1e-2,1e-3];
    fs = []
    for x0=xs
        for t0=ts
            disp(x0)
            [prob,data,opts,h] = integrator(false,N_sim, x0);
            x_curr = [x0;t0];
            x_res = x_curr;
            lambda_curr = 0;
            cost = 0;
            for step=1:N_sim
                prob.w.x(0,0,data.n_s).init = x_curr;
                prob.w.x(0,0,data.n_s).lb = x_curr;
                prob.w.x(0,0,data.n_s).ub = x_curr;
                prob.w.lambda(0,0,data.n_s).init = lambda_curr;
                prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
                prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
                success = homotopy(prob, comp_tol, comp_tol);
                if ~success
                    disp(['Failure to converge at step=' num2str(step)])
                else
                    disp(['step=' num2str(step)])
                end
                x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
                lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
                x_sim = prob.w.x(1,:,data.n_s).res;
                x_res = [x_res,x_sim];
                %disp(prob.g.objective(1).res)
                %cost = cost + prob.f_result;
                cost = cost + prob.g.objective(1).res;
            end
            x_curr;
            x_term = [x_term,x_curr];
            f = cost + (x_curr(1)-3.5)^2;
            disp(f)
            fs = [fs,f];
        end
    end
    plot(xs, fs, "DisplayName", ['$\sigma =' num2str(comp_tol) '$'])
end
