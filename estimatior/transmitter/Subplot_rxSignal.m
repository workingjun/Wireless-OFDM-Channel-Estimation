function []=Subplot_rxSignal(params)
    GP = 16;
    NN = 1:GP;
    
    fig = figure(101);
    fig.Position = [620, 0, 500, 600];
    sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);

    subplot(3, 1, 1);
    stem(NN, params.e_rx); 
    xlim([0 17]);
    title("R.V epcilon");
   
    subplot(3, 1, 2); 
    stem(NN, params.e_rx_pdf);
    xlim([0 17]);
    title('f_p');

    subplot(3, 1, 3);         
    stem(NN, params.e_rx_multi_pdf);
    xlim([0 17]);
    title(['LR of Joint $f_e$ / $\hat{L}$: '  num2str(params.L_sol_erx)], 'Interpreter', 'latex');
   
    pause;
end