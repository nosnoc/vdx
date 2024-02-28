%%

stage_counts = [10, 50, 100];%, 500, 1000, 5000, 10000];
n_s_list = [1,2,3,4];
N_fe = 2;
T = (11/12)*pi + sqrt(3);

% do no FESD runs
no_fesd_all_hs = {};
no_fesd_all_errors = {};
for n_s=n_s_list
    no_fesd_errors = [];
    no_fesd_hs = [];
    for N_sim=stage_counts
        [prob, data, opts, h] = integrator(T, N_sim, N_fe, false, n_s);
        orig_init = prob.w.init;
        %% S I M U L A T E
        x_curr = data.x0;
        lambda_curr = 0;
        for step=1:N_sim
            prob.w.init = orig_init;
            prob.w.x(0,0,data.n_s).init = x_curr;
            prob.w.x(0,0,data.n_s).lb = x_curr;
            prob.w.x(0,0,data.n_s).ub = x_curr;
            prob.w.lambda(0,0,data.n_s).init = lambda_curr;
            prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
            prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
            [success,stats] = homotopy(prob, 1, 1e-12);
            disp(['step=' num2str(step)])
            x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
            lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
        end
        error = norm(x_curr - [-1;0]);
        no_fesd_hs = [no_fesd_hs;h];
        no_fesd_errors = [no_fesd_errors;error];
    end
    no_fesd_all_hs{n_s} = no_fesd_hs;
    no_fesd_all_errors{n_s} = no_fesd_errors;
end

% do FESD runs
fesd_all_hs = {};
fesd_all_errors = {};
for n_s=n_s_list
    fesd_errors = [];
    fesd_hs = [];
    for N_sim=stage_counts
        [prob, data, opts, h] = integrator(T, N_sim, N_fe, true, n_s);
        orig_init = prob.w.init;
        %% S I M U L A T E
        x_curr = data.x0;
        lambda_curr = 0;
        for step=1:N_sim
            prob.w.init = orig_init;
            prob.w.x(0,0,data.n_s).init = x_curr;
            prob.w.x(0,0,data.n_s).lb = x_curr;
            prob.w.x(0,0,data.n_s).ub = x_curr;
            prob.w.lambda(0,0,data.n_s).init = lambda_curr;
            prob.w.lambda(0,0,data.n_s).lb = lambda_curr;
            prob.w.lambda(0,0,data.n_s).ub = lambda_curr;
            [success,stats] = homotopy(prob, 1, 1e-12);
            disp(['step=' num2str(step)])
            x_curr = prob.w.x(1,data.N_fe,data.n_s).res;
            lambda_curr = prob.w.lambda(1,data.N_fe,data.n_s).res;
        end
        error = norm(x_curr - [-1;0]);
        fesd_hs = [fesd_hs;h];
        fesd_errors = [fesd_errors;error];
    end
    fesd_all_hs{n_s} = fesd_hs;
    fesd_all_errors{n_s} = fesd_errors;
end

%%plot
figure
hold on;

for n_s=n_s_list
    plot(no_fesd_all_hs{n_s},no_fesd_all_errors{n_s}, 'DisplayName', sprintf('Standard $n_s=%d$', n_s));
    plot(fesd_all_hs{n_s},fesd_all_errors{n_s}, 'DisplayName', sprintf('FESD $n_s=%d$', n_s));
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
