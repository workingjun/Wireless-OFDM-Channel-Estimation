function []= Subplot_method2_upgrade_LLR(params, ylim_flag)
    fig = figure(1);
    fig.Position = [300, 0, 900, 600];
    sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);

    % 첫 번째 행
    subplot(2,3,1);
    stem(1:16, params.p_rx_sum_pdf1);
    title(['figure101 $\hat{L}$: ' num2str(params.p_sol1)], 'Interpreter', 'latex');
    if ylim_flag && sum(params.p_rx_sum_pdf1) > 0
       ylim([0 30]);
    end
    grid on;
    
    subplot(2,3,2);
    stem(1:16, params.e_rx_sum_pdf1);
    title(['figure102 $\hat{L}$: ' num2str(params.e_sol1)], 'Interpreter', 'latex');
    ylim([0 30]);

    grid on;
    
    subplot(2,3,3);
    stem(1:16, params.pe_rx_sum_pdf1);
    title(['figure103 $\hat{L}$: ' num2str(params.pe_sol1)], 'Interpreter', 'latex');
    if ylim_flag && sum(params.pe_rx_sum_pdf1) > 0
       ylim([0 30]);
    end
    grid on;
    
    % 두 번째 행
    subplot(2,3,4);
    stem(1:16, params.p_rx_sum_pdf2);
    title(['figure1001 $\hat{L}$: ' num2str(params.p_sol2)], 'Interpreter', 'latex');
    % ylim([params.p_rx_sum_pdf2(4)-100 params.p_rx_sum_pdf2(4)+100]);
    
    grid on;
    
    subplot(2,3,5);
    stem(1:16, params.e_rx_sum_pdf2);
    title(['figure1002 $\hat{L}$: ' num2str(params.e_sol2)], 'Interpreter', 'latex');
    if ylim_flag
       ylim([0 max(params.e_rx_sum_pdf2)+10]);
    end
    % ylim([-10 0]);
    
    grid on;
    
    subplot(2,3,6);
    stem(1:16, params.pe_rx_sum_pdf2);
    title(['figure1003 $\hat{L}$: ' num2str(params.pe_sol2)], 'Interpreter', 'latex');
    
    % ylim([max(params.pe_rx_sum_pdf2)-30 max(params.pe_rx_sum_pdf2)+10]);
    
    grid on;
    
    pause; 

end