function []=PlotBER_SER(params)
    
    fig = figure(2);
    fig.Position = [0, 0, 800, 600];

    semilogy(params.SNR_dB, params.Sim_BER, 'b+-');%% Simulation !!
    hold on;
    semilogy(params.SNR_dB, params.Sim_SER, 'ro-');%% Simulation !!
    semilogy(params.SNR_dB, params.Sim_BER2, 'bs-');%% Simulation !!
    semilogy(params.SNR_dB, params.Sim_SER2, 'rh-');%% Simulation !!
    semilogy(params.SNR_dB, Analysis_SER,'kx-');%%  성능분석 결과
    semilogy(params.SNR_dB, Analysis_BER,'k*-');%%  성능분석 결과
    legend('BER1, Simulation', 'SER1, Simulation' , 'BER2, Simulation', 'SER2, Simulation' , 'Fading BER, Analysis','AWGN BER, Analysis', 'Location','northeast');
    grid on;
    xlabel('SNR (dB)', 'Fontsize', 12);
    ylabel('Error Rate', 'Fontsize', 12);
    title('MPSK (Fading Channel) [by KKB]', 'Fontsize', 14);
    xlim([min(params.SNR_dB) max(params.SNR_dB)]);
    ylim([10^-7 1]);

end