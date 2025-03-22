function []=Subplot_method2_upgrade(params, flags)  
    if flags
        fig = figure(102);
        fig.Position = [200, 0, 1280, 800];

        for u=1:params.GP
            subplot(4, 8, 2*u-1); 
            stem(1:16, params.SIM_channel_part(u, :));
            grid on;
            xlim([0 17]);
            ylim([0 5]);
            title(['channel / u=' num2str(u)], 'Interpreter', 'latex');
    
            subplot(4, 8, 2*u); 
            stem(1:16, params.SIM_noise_part(u, :));
            grid on;
            xlim([0 17]);
            ylim([-10 10]);
            title(['noise / u=' num2str(u)], 'Interpreter', 'latex');
        end
       
    else
        fig = figure(101);
        fig.Position = [200, 100, 1000, 600];
        sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);
    
        subplot(2, 3, 1); 
        stem(1:16, params.ABdiffsq);
        grid on;
        xlim([0 17]);
        title('Random Variable / 1/N ', 'Interpreter', 'latex');
        
        subplot(2, 3, 2); 
        stem(1:16, params.ABdiffsq_ch);
        grid on;
        xlim([0 17]);
        title('Random Variable / not 1/N ', 'Interpreter', 'latex');
    
        subplot(2, 3, 3); 
        stem(1:16, params.e_sum_pdf);
        grid on;
        xlim([0 17]);
        ylim([-5 10]);
        title(['Method2 / $\hat{L}$: '  num2str(params.my_method2_sol)], 'Interpreter', 'latex');
    
        subplot(2, 3, 4);
        stem(1:16, params.SIM_channel); 
        grid on;
        xlim([0 17]);
        ylim([-5 10]);
        title(['channel part / $\hat{L}$: ' num2str(params.confirm_sol)], 'Interpreter', 'latex');
       
        subplot(2, 3, 5);
        stem(1:16, params.SIM_noise);
        grid on;
        xlim([0 17]);
        ylim([-5 10]);
        title(['noise part / $\hat{L}$: ' num2str(params.confirm2_sol)], 'Interpreter', 'latex');
    
        subplot(2, 3, 6);         
        stem(1:16, params.SIM_upgrade);
        grid on;
        xlim([0 17]);
        ylim([-5 20]);
        title(['Method2 upgrade / $\hat{L}$: '  num2str(params.upgrade_sol)], 'Interpreter', 'latex');
    end

    pause;
end