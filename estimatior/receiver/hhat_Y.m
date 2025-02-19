function [Y_pdf, Y_multi_pdf, Y_sum_pdf, L_sol_y1, L_sol_y2] = hhat_Y(params)
    N_OFDM_symbols = 100;
    N = params.N;
    L = params.L;
    
    mean_z = params.JJ(L+1); 
    calc_z = params.JJ(L+1)^2;
    sigma_z = calc_z / N_OFDM_symbols;

    Y_pdf = (2 * pi * sigma_z)^(-0.5) * exp(-0.5 * ((params.Y_u - mean_z).^2 / sigma_z));
    Y_log_pdf = -0.5 * ((params.Y_u - mean_z).^2 / sigma_z) - log(sqrt(2 * pi * sigma_z));

    Y_multi_pdf = cumprod(Y_pdf, "reverse");
    Y_sum_pdf = cumsum(Y_log_pdf, "reverse");
    
    scale_factors = 1 ./ (N + 1 - (1:N));
    Y_multi_pdf(1:N) = Y_multi_pdf(1:N) .^ scale_factors;
    Y_sum_pdf(1:N) = Y_sum_pdf(1:N) .* scale_factors;
    
    [~, L_sol_y1] = max(Y_multi_pdf(1:16));
    [~, L_sol_y2] = max(Y_sum_pdf(1:16));

    L_sol_y1 = L_sol_y1 - 1;
    L_sol_y2 = L_sol_y2 - 1;
end
