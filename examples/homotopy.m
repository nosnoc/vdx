function success = homotopy(prob,slope)
    if ~exist('slope')
        slope = 0.1;
    end
    sigma_k = 1;
    all_stats = [];
    comp_tol = 1e-7;
    while sigma_k >= comp_tol
        prob.p.sigma(1).init = sigma_k;
        stats = prob.solve();
        prob.w.init = prob.w.res;
        all_stats = [all_stats, stats];
        comp_res = full(prob.comp_res_fun(prob.w.res,prob.p.init));
        fprintf(['sigma_k=' num2str(sigma_k) ' comp_res=' num2str(comp_res) '\n']);
        sigma_k = slope*sigma_k;
        if comp_res < comp_tol
            break
        end
    end

    success = stats.success;
    if ~success
        1+1
    end
end
