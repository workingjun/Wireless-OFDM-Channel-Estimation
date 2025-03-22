function [Rx_Bits] = Rx_Bits_gen(M, rx_signal)
%UNTITLED2 이 함수의 요약 설명 위치
%   자세한 설명 위치

        if M == 2 
            if real(rx_signal) > 0 %% Bits Estimation :: 수신 비트 검출(결정) !!
                Rx_Bits(1) = 1; 
            else
                Rx_Bits(1) = 0;
            end
        elseif M == 4
            if real(rx_signal) > 0 %% Bits Estimation :: 수신 비트 검출(결정) !!
                Rx_Bits(1) = 1; 
            else
                Rx_Bits(1) = 0;
            end
            if imag(rx_signal) > 0 %% Bits Estimation :: 수신 비트 검출(결정) !!
                Rx_Bits(2) = 1; 
            else
                Rx_Bits(2) = 0;
            end
        elseif M == 8
            if real(rx_signal) > 0 %% Bits Estimation :: 수신 비트 검출(결정) !!
                Rx_Bits(2) = 0; 
            else
                Rx_Bits(2) = 1;
            end
            if imag(rx_signal) > 0 %% Bits Estimation :: 수신 비트 검출(결정) !!
                Rx_Bits(1) = 0; 
            else
                Rx_Bits(1) = 1;
            end
            if abs(angle(rx_signal))>pi/4 && abs(angle(rx_signal))<pi*3/4
                Rx_Bits(3) = 1; 
            else
                Rx_Bits(3) = 0; 
            end
            %%% 8PSK에 대해 여러분이 코딩을 해야됨!
            %%% Out_Bits(1) = ???
            %%% Out_Bits(2) = ???
            %%% Out_Bits(3) = ???
        end



end

