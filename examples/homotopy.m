function success = homotopy(prob,slope)
    if ~exist('slope')
        slope = 0.1;
    end
    sigma_k = 1;
    all_stats = [];
    while sigma_k >= 1e-12
        prob.p.sigma(1).init = sigma_k;
        stats = prob.solve();
        sigma_k = slope*sigma_k;
        prob.w.init = prob.w.res;
        all_stats = [all_stats, stats];
    end

    success = stats.success;
    if ~success
        1+1
    end
end
