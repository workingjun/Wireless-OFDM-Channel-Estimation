function []=Ber_Performance(params)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN 성능 분석 결과!
    Analysis_BER = zeros(1,length(params.SNR_dB));
    Analysis_SER = zeros(1,length(params.SNR_dB));
    if params.M == 2 
        Analysis_BER = 1/2*erfc( sqrt(params.SNR) );
        %Analysis_SER = 1/2*erfc( sqrt(params.SNR) );
        Analysis_SER = 1/2*(1 - sqrt( params.SNR./(params.SNR+1)));
    elseif params.M == 4 
        Analysis_BER = 1/2*erfc( sqrt(params.SNR/2) );%% 왜 SNR/2 인가?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SER_app = 2*Q( sqrt(2*NSR)*sin(pi/params.M) );
        %Analysis_SER = 1*erfc( sqrt(params.SNR)*sin(pi/params.M) );%%%%%%%%%%%%%%%%%%%%%%%%%% 근사화된 식!!!
        Analysis_SER = 1/2.*(1 - sqrt( (params.SNR/2)./((params.SNR/2)+1)));
    elseif params.M == 8 
        Analysis_SER = 1*erfc( sqrt(params.SNR)*sin(pi/params.params.M) );%%%%%%%%%%%%%%%%%%%%%%%%%% 근사화된 식!!!
        Analysis_BER = (1/log2(params.M))*Analysis_SER;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 근사화된 식!!!
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 성능 분석 결과!
    
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