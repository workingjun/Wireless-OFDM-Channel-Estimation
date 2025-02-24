function [ABdiffsq, ABdiffsq_ch, mean_p, ...
          p_rx_logpdf, p_rx_sum_pdf, e_rx_sum_pdf, ...
          SIM_upgrade, SIM_channel, SIM_noise, ...
          SIM_channel_part, SIM_noise_part, ...
          upgrade_sol, confirm_sol, confirm2_sol] = method2_upgraded(params)

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
     
    mean_p = (params.c_hat*(1-params.rho(4,16:-1:1))+params.J)./(2*(GP+1-(1:GP))); %%평균
    sigma_p = mean_p.^2 / N_OFDM_symbols; %%분산
    p_rx_logpdf = -0.5 * (abs(ABdiffsq - mean_p).^2 ./ sigma_p) - log(sqrt(2 * pi * sigma_p));

    % mean_p = (params.c_hat*(1-params.rho(4,16:-1:1))+params.J)/2; %%평균
    % sigma_p = mean_p.^2 / N_OFDM_symbols; %%분산
    % p_rx_logpdf = -0.5 * (abs(ABdiffsq_ch - mean_p).^2 ./ sigma_p) - log(sqrt(2 * pi * sigma_p));
    % 
    p_rx_sum_pdf = cumsum(p_rx_logpdf, 'reverse');
    p_rx_sum_pdf = p_rx_sum_pdf .* (1 ./ (GP + 1 - (1:GP)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % mean_e = params.J(4) ./ (GP-(1:GP)+1); %%평균
    % sigma_e = mean_e.^2 / N_OFDM_symbols; %%분산
    % 
    % e_rx_logpdf = -0.5 * (abs(ABdiffsq - mean_e).^2 ./ sigma_e) - log(sqrt(2 * pi * sigma_e));

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

    [~, upgrade_sol] = max(p_rx_sum_pdf);
    [~, confirm_sol] = max(SIM_channel);
    [~, confirm2_sol] = max(SIM_noise);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



