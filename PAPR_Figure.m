function [out] = PAPR_Figure(tx_signal, GP, N)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
%
% h :: 시간축 채널 값
% H :: 주파수 축 채널 값

    No_OFDM_symbol = (length(tx_signal))/(GP+N);

    t_index=1:N;
    for ss=1:No_OFDM_symbol
        ST = (ss-1)*(GP+N) + GP+1;
        To = (ss-1)*(GP+N) + (GP+N); 

        t_signal = tx_signal(ST:To);
        
        Pw = abs(t_signal).^2;
        PAPR = max(Pw) / (sum(Pw)/N);
        
        [max(Pw) (sum(Pw)/N) PAPR]

        figure(2000);
        semilogy(t_index, abs(t_signal).^2);
        hold on;
        grid on;
        xlabel('Time (GP+N)', 'Fontsize', 12);
        ylabel('|x(n)|^2 (dB)', 'Fontsize', 12);
        title('PAPR [by KKB]', 'Fontsize', 14);
        xlim([min(t_index) max(t_index)]);
        ylim([10^-2 10]);
        pause();


    end
    out = 1;


end

