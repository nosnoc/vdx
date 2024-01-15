function homotopy(prob)
    sigma_k = 1;
    while sigma_k >= 1e-7
        prob.p.sigma(1).init = sigma_k;
        prob.solve();
        sigma_k = 0.1*sigma_k;
        prob.w.init = prob.w.res;
    end
end
