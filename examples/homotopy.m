function success = homotopy(prob,slope)
    if ~exist('slope')
        slope = 0.1;
    end
    sigma_k = 1;
    while sigma_k >= 1e-12
        prob.p.sigma(1).init = sigma_k;
        stats = prob.solve();
        sigma_k = slope*sigma_k;
        prob.w.init = prob.w.res;
    end

    success = stats.success;
end
