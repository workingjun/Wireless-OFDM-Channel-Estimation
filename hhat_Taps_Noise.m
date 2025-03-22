function [p_pdf, z_pdf, p_Log_pdf, z_Log_pdf, L_sol_p, L_sol_z] = hhat_Taps_Noise(params)
    N_OFDM_symbols = 100;
    N = params.N;
    L = params.L;
    
    mean_z = params.JJ(L+1); %%평균
    calc_z = params.JJ(L+1)^2/(N-L); %%분산의 분자 계산 
    sigma_z = calc_z / N_OFDM_symbols; %%분산   
    
    h_sq = sum(params.Y_u(1:N)) - N*params.JJ(L+1); 

    mean_p = h_sq + L * params.JJ(L+1); %%평균
    calc_p = L * params.JJ(L+1)^2 + 2 * h_sq * params.JJ(L+1); %%분산의 분자 계산
    sigma_p = calc_p / N_OFDM_symbols; %%분산

    p_Log_pdf = -0.5 * ((params.P(1:N-1) - mean_p).^2 / sigma_p) - log(sqrt(2 * pi * sigma_p));
    z_Log_pdf = -0.5 * ((params.Z(1:N-1) - mean_z).^2 / sigma_z) - log(sqrt(2 * pi * sigma_z));

    p_pdf = (2 * pi * sigma_p)^(-0.5) * exp(-0.5 * ((params.P(1:N-1) - mean_p).^2 / sigma_p));
    z_pdf = (2 * pi * sigma_z)^(-0.5) * exp(-0.5 * ((params.Z(1:N-1) - mean_z).^2 / sigma_z));

    [~, L_sol_p] = max(p_pdf(1:16));
    [~, L_sol_z] = max(z_pdf(1:16));
end
