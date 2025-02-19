function [] = Subplot_hhat(params)
    N = 64;
    NN = 1:N;
        
    fig = figure(101);
    fig.Position = [0, 0, 900, 900];
    sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);

    subplot(3, 3, 1);
    stem(NN, params.Y_pdf); 
    xlim([0 65]);
    title("f_Y");

    subplot(3, 3, 4); 
    stem(NN, params.Y_multi_pdf);
    xlim([0 65]);
    title(['LR of Joint $f_Y$ / $\hat{L}$: '  num2str(params.L_sol_y1)], 'Interpreter', 'latex');

    subplot(3, 3, 7); 
    stem(NN, params.Y_sum_pdf);
    xlim([0 65]);
    title(['LLR of Joint $f_Y$ / $\hat{L}$: '  num2str(params.L_sol_y2)], 'Interpreter', 'latex');
    ylim([-1 15]);

    subplot(3, 3, 2);
    stem(NN, params.e_pdf);
    xlim([0 65]);
    title("f_e");

    subplot(3, 3, 5);
    stem(NN, params.e_multi_pdf);
    xlim([0 65]);
    title(['LR of Joint $f_e$ / $\hat{L}$: '  num2str(params.L_sol_e1)], 'Interpreter', 'latex');

    subplot(3, 3, 8);
    stem(NN, params.e_sum_pdf);
    xlim([0 65]);
    title(['LLR of Joint $f_e$ / $\hat{L}$: '  num2str(params.L_sol_e2)], 'Interpreter', 'latex');
    ylim([-1 15]);

    subplot(3, 3, 3);
    stem(1:N-1, params.p_pdf);
    xlim([0 64]);
    title('f_p');    

    subplot(3, 3, 6);
    stem(NN, params.Pe_multi_pdf);
    xlim([0 65]);
    title(['LR of Joint $f_{Pe}$ / $\hat{L}$: '  num2str(params.L_sol_pe1)], 'Interpreter', 'latex');   

    subplot(3, 3, 9);
    stem(NN, params.Pe_sum_pdf);
    xlim([0 65]);
    title(['LLR of Joint $f_{Pe}$ / $\hat{L}$: '  num2str(params.L_sol_pe2)], 'Interpreter', 'latex');   
    ylim([-1 15]);

    pause;  
end
            