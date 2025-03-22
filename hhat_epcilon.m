function [e_pdf, e_multi_pdf, e_sum_pdf, L_sol_e1, L_sol_e2] = hhat_epcilon(params)
    
    N = params.N;
    L = params.L;
    e_u = zeros(1, N);
    N_OFDM_symbols = 100;
    
    e_u(1:N-1) = params.JJ(1:N-1) - params.JJ(2:N) .* (1 - 1 ./ (N - (1:N-1) + 1));
    e_u(N) = params.JJ(N);
    
    mean_e = params.JJ(L+1) ./ (N - (1:N) + 1);
    calc_e = params.JJ(L+1)^2 ./ (N - (1:N) + 1).^2;
    sigma_e = calc_e / N_OFDM_symbols;
    
    e_pdf = (2 * pi * sigma_e).^(-0.5) .* exp(-0.5 * ((e_u - mean_e).^2 ./ sigma_e));
    e_log_pdf = -0.5 * ((e_u - mean_e).^2 ./ sigma_e) - log(sqrt(2 * pi * sigma_e));
    
    e_multi_pdf = cumprod(e_pdf, "reverse");
    e_sum_pdf = cumsum(e_log_pdf, "reverse");
    
    scale_factors = 1 ./ (N - (1:N) + 1);
    e_multi_pdf = e_multi_pdf .^ scale_factors;
    e_sum_pdf = e_sum_pdf .* scale_factors;
    
    [~, L_sol_e1] = max(e_multi_pdf(1:16));
    [~, L_sol_e2] = max(e_sum_pdf(1:16));
    
    L_sol_e1 = L_sol_e1 - 1;
    L_sol_e2 = L_sol_e2 - 1;
end
