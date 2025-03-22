function [e_sol1, e_sol2, new_sol2] = method2_upgraded(params)
    GP = params.GP;
    L = params.L;
    N = params.N;
    r = params.rx_signal;
    N_OFDM_symbols = length(r)/(GP+N);
    
    ABdiffsq = params.ABdiffsq;
    ABdiffsq_ratio = params.ABdiffsq_ratio;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % mean_p = zeros(1, GP);
    % for kk=1:GP
    %     mean_p(kk) = (params.c_hat*(1-params.rho(kk, GP-kk+1)));
    % end
    % 
    % sigma_p = mean_p.^2 / N_OFDM_symbols; 
    % p_rx_logpdf = -0.5 * (abs(ABdiffsq_ch - mean_p).^2 ./ sigma_p) - log(sqrt(2 * pi * sigma_p));
    % 
    % p_rx_sum_pdf1 = cumsum(p_rx_logpdf);
    % p_rx_sum_pdf1 = p_rx_sum_pdf1./(1:GP);
    % 
    % p_rx_sum_pdf2 = cumsum(p_rx_logpdf);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ep_ratio = ABdiffsq_ratio;
    % ep = ABdiffsq;

    ep_ratio_LLRs = zeros(GP, GP);
    ep_ratio_PDFs = zeros(GP, GP);
    ep_ratio_LLR_sum2 = zeros(1, GP);
    ep_ratio_PDF_prod2 = ones(1, GP);

    for u=1:GP
        for k=u:GP
            m_e = params.J(u)/(GP-k+1);
            var_e = m_e^2/N_OFDM_symbols;
            ep_ratio_LLRs(u, k) = -abs(ep_ratio(k)-m_e)^2/(2*var_e)-log(2*pi*var_e)/2;
            ep_ratio_PDFs(u, k) = exp(ep_ratio_LLRs(u, k)); 
        end
        ep_ratio_LLR_sum2(u) = sum(ep_ratio_LLRs(u, u:GP))/(GP-u+1);
        ep_ratio_PDF_prod2(u) = prod(ep_ratio_PDFs(u, u:GP))^(1/(GP-u+1));
    end

    % for uu=1:GP
    %     stem(1:16, ep_LLRs(uu,:)); grid on;
    %     pause()
    % end
    % figure(11)
    % stem(1:16, ep_ratio_LLR_sum); grid on;
    % % ylim([-10 100]);
    % pause()
    % figure(12)
    % stem(1:16, ep_ratio_LLR_sum2); grid on;
    % % ylim([-10 20]);
    % pause()
    % figure(13);
    % stem(1:16, ep_ratio_PDF_prod2); grid on;
    % pause()


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ep = ABdiffsq;
    
    ep_LLRs = zeros(1, GP);
    ep_PDFs = zeros(1, GP);
    ep_LLR_sum = zeros(1, GP);
    ep_LLR_sum2 = zeros(1, GP);
    ep_PDF_prod = ones(1, GP);
    ep_PDF_prod2 = ones(1, GP);
    
    %%% uu 가 고정되면, 'GP-uu+1'개의 LLR 계산하는 동안에 J(uu)만을 사용!!
    for u=1:GP %% For the given uu, J(uu) is fixed !!
        m_e = params.J(u);
        var_e = m_e^2/N_OFDM_symbols;

        for k=u:GP
            % m_e = params.J(k);
            % var_e = m_e^2/N_OFDM_symbols;
            ep_LLRs(u, k) = -abs(ep(k)-m_e)^2/(2*var_e)-log(2*pi*var_e)/2;
            ep_PDFs(u, k) = exp(ep_LLRs(u, k)); 
        end

        ep_LLR_sum(u) = sum(ep_LLRs(u, u:GP));
        ep_PDF_prod(u) = prod(ep_PDFs(u, u:GP));
        ep_LLR_sum2(u) = ep_LLR_sum(u)/(GP-u+1);
        ep_PDF_prod2(u) = ep_PDF_prod(u)^(1/(GP-u+1));

    end
    
    % figure(10)
    % for uu=1:GP
    %     stem(1:16, ep_LLRs(uu,:)); grid on;
    %     pause()
    % end
    % figure(11)
    % stem(1:16, ep_LLR_sum); grid on;
    % % ylim([-10 100]);
    % pause()
    % figure(12)
    % stem(1:16, ep_LLR_sum2); grid on;
    % % ylim([-10 20]);
    % pause()
    % figure(13);
    % stem(1:16, ep_PDF_prod); grid on;
    % pause()
    % figure(14);
    % stem(1:16, ep_PDF_prod2); grid on;
    % pause()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % fig = figure(1);
    % fig.Position = [300, 0, 600, 800];
    % sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);
    % 
    % subplot(3, 2, 1)
    % stem(1:16, ep_ratio_LLR_sum2);
    % ylim([0 20])
    % grid on;
    % 
    % subplot(3, 2, 2)
    % stem(1:16, ep_ratio_PDF_prod2);
    % grid on;
    % 
    % subplot(3, 2, 3)
    % stem(1:16, ep_LLR_sum);
    % ylim([-10 100]);
    % grid on;
    % 
    % subplot(3, 2, 4)
    % stem(1:16, ep_PDF_prod);
    % grid on;
    % 
    % subplot(3, 2, 5)
    % stem(1:16, ep_LLR_sum2);
    % ylim([-2 10]);
    % grid on;
    % 
    % subplot(3, 2, 6)
    % stem(1:16, ep_PDF_prod2);
    % grid on;
    % 
    % pause;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ep = ABdiffsq;

    ep_channel_LLRs = zeros(GP, GP);
    ep_noise_LLRs = zeros(GP, GP);
    ep_ch_LLR_sum = zeros(1, GP);
    ep_no_LLR_sum = zeros(1, GP);
    ep_total_sum = zeros(1, GP);
  
    for u=1:GP
        m_p = params.c_hat*(1-params.rho(u, GP-u+1));
        var_p = m_p^2/N_OFDM_symbols;
            
        m_e = params.J(u);
        var_e = m_e^2/N_OFDM_symbols;

        for k=1:u
            ep_channel_LLRs(u, k) = -abs(ep(k)-m_p)^2/(2*var_p)-log(2*pi*var_p)/2;
        end
        for k=u:GP    
            ep_noise_LLRs(u, k) = -abs(ep(k)-m_e)^2/(2*var_e)-log(2*pi*var_e)/2;
        end

        ep_ch_LLR_sum(u) = sum(ep_channel_LLRs(u, 1:u));
        ep_no_LLR_sum(u) = sum(ep_noise_LLRs(u, u:GP));
        ep_total_sum(u) = ep_ch_LLR_sum(u) + ep_no_LLR_sum(u);
    end
    
    % figure(10)
    % for uu=1:GP
    %     disp(ep_channel_LLRs(uu,:));
    %     stem(1:16, ep_channel_LLRs(uu,:)); grid on;
    %     ylim([-10 10]);
    %     pause()
    % end
    % 
    % figure(11)
    % for uu=1:GP
    %     stem(1:16, ep_noise_LLRs(uu,:)); grid on;
    %     pause()
    % end

    % figure(12)
    % stem(1:16, ep_ch_LLR_sum); grid on;
    % % ylim([-500 10]);
    % pause()
    % 
    % figure(13)
    % stem(1:16, ep_no_LLR_sum); grid on;
    % % ylim([-10 100]);
    % pause()

    figure(14)
    stem(1:16, ep_total_sum); grid on;
    % ylim([-1000 0]);
    [~, ppp] = max(ep_total_sum);
    title(['$\hat{L}$: ' num2str(ppp)], 'Interpreter', 'latex');
    pause()
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, e_sol1] = max(ep_LLR_sum);
    [~, new_sol2] = max(ep_ratio_LLR_sum2);
    [~, e_sol2] = max(ep_LLR_sum2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



