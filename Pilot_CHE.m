function [hat_H] = Pilot_CHE(No_Pilot_symbols, N, GP, rx_signal, Tx_Symbols)
%UNTITLED �� �Լ��� ��� ���� ��ġ

    hat_H = zeros(1,N);
    for tt=1:No_Pilot_symbols
        ST = (tt-1)*(N+GP) + GP+1; 
        To = (tt-1)*(N+GP) + GP+N; 
        Rx_Signal_T = rx_signal(ST:To);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ð��࿡�� GP�� ������ OFDM �ɺ�!!
        %
        temp = fft(Rx_Signal_T)/sqrt(N);
        Rx_Signal_F = fftshift(temp);

        ST = (tt-1)*N + 1; 
        To = (tt-1)*N + N; 
        Pilots = Tx_Symbols(ST:To);

        hat_H = hat_H + Rx_Signal_F./Pilots;
    end
    hat_H = hat_H/No_Pilot_symbols; 

end

