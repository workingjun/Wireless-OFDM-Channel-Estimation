function [p_rx_logpdf1, p_rx_logpdf2, p_rx_sum_pdf1, p_rx_sum_pdf2, ...
          e_rx_logpdf1, e_rx_logpdf2, e_rx_sum_pdf1, e_rx_sum_pdf2, ...
          pe_rx_sum_pdf1, pe_rx_sum_pdf2, ...
          p_sol1, p_sol2, e_sol1, e_sol2, pe_sol1, pe_sol2] = method2_upgraded(params)

    GP = params.GP;
    L = params.L;
    N = params.N;
    r = params.rx_signal;
    N_OFDM_symbols = length(r)/(GP+N);

    SIM_upgrade = zeros(1, GP-1);
    e_rx_sum_pdf = zeros(GP, GP);
    p_rx_sum_pdf = zeros(GP, GP);
    SIM_noise = zeros(1, GP);
    SIM_channel = zeros(1, GP);
    SIM_noise_part = zeros(GP, GP);
    SIM_channel_part = zeros(GP, GP);
    
    ABdiffsq = params.ABdiffsq;
    ABdiffsq_ch = params.ABdiffsq_ch;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mean_p = (params.c_hat*(1-arrayfun(@(idx) params.rho(idx, GP-idx+1), 1:GP))+params.J)./(2*(GP+1-(1:GP)));
    sigma_p = mean_p.^2 / N_OFDM_symbols; 
    p_rx_logpdf1 = -0.5 * (abs(ABdiffsq - mean_p).^2 ./ sigma_p) - log(sqrt(2 * pi * sigma_p));
    
    mean_p = (params.c_hat*(1-arrayfun(@(idx) params.rho(idx, GP-idx+1), 1:GP))+params.J)/2;
    sigma_p = mean_p.^2 / N_OFDM_symbols; 
    p_rx_logpdf2 = -0.5 * (abs(ABdiffsq_ch - mean_p).^2 ./ sigma_p) - log(sqrt(2 * pi * sigma_p));

    p_rx_sum_pdf1 = cumsum(p_rx_logpdf1, "reverse");
    p_rx_sum_pdf1 = p_rx_sum_pdf1./(1:GP);

    p_rx_sum_pdf2 = cumsum(p_rx_logpdf2);
    p_rx_sum_pdf2 = p_rx_sum_pdf2./(1:GP);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mean_e = params.J ./ (GP-(1:GP)+1);
    sigma_e = mean_e.^2 / N_OFDM_symbols; 
    e_rx_logpdf1 = -0.5 * (abs(ABdiffsq - mean_e).^2 ./ sigma_e) - log(sqrt(2 * pi * sigma_e));

    mean_e = params.J;
    sigma_e = mean_e.^2 / N_OFDM_symbols; 
    e_rx_logpdf2 = -0.5 * (abs(ABdiffsq_ch - mean_e).^2 ./ sigma_e) - log(sqrt(2 * pi * sigma_e));
    
    e_rx_sum_pdf1 = cumsum(e_rx_logpdf1, "reverse");
    e_rx_sum_pdf1 = e_rx_sum_pdf1./(GP+1-(1:GP));

    e_rx_sum_pdf2 = cumsum(e_rx_logpdf2, "reverse");
    e_rx_sum_pdf2 = e_rx_sum_pdf2./(GP+1-(1:GP));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pe_rx_sum_pdf1 = p_rx_sum_pdf1 + e_rx_sum_pdf1;
    pe_rx_sum_pdf2 = p_rx_sum_pdf2 + e_rx_sum_pdf2;
    
    [~, p_sol1] = max(p_rx_sum_pdf1);
    [~, p_sol2] = max(p_rx_sum_pdf2);
    [~, e_sol1] = max(e_rx_sum_pdf1);
    [~, e_sol2] = max(e_rx_sum_pdf2);
    [~, pe_sol1] = max(pe_rx_sum_pdf1);
    [~, pe_sol2] = max(pe_rx_sum_pdf2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 이전에 했던 방식; 노이즈 파트가 다름
    % p_rx_sum_pdf = zeros(1, GP-1); e_rx_sum_pdf = zeros(1, GP-1); 
    % for u = 1:GP-1
    %     p_rx_sum_pdf(u) = sum(p_rx_logpdf(1:u))/u;
    %     e_rx_sum_pdf(u) = sum(e_rx_logpdf(u+1:GP))/(GP-u);
    %     SIM_upgrade(u) =  p_rx_sum_pdf(u) + e_rx_sum_pdf(u);
    % end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for u = 1:GP
    %     SIM_channel_part(u, 1:u) = p_rx_logpdf(1:u);
    %     SIM_noise_part(u, u:GP) = e_rx_logpdf(u:GP);
    %     SIM_channel(u) = sum(SIM_channel_part(u, :))/(u);
    %     SIM_noise(u) = sum(SIM_noise_part(u, :))/(GP-u+1);
    % 
    %     % SIM_channel_part(u, u) = p_rx_logpdf(u);
    %     % SIM_noise_part(u, u:GP) = e_rx_logpdf(u:GP);
    %     % SIM_channel(u) = sum(SIM_channel_part(u, :));
    %     % SIM_noise(u) = sum(SIM_noise_part(u, :))/(GP-u+1);
    % 
    %     SIM_upgrade(u) = SIM_channel(u) + SIM_noise(u);
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



