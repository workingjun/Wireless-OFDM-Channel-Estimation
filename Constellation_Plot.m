function [out] = Constellation_Plot(M, Tx_Bits, rx_signal)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
%   [Tx_Bits(1),Tx_Bits(2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [0,1] | [1,1]
% --------------------------
%   [0,0] | [1,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if M == 4
%     if Tx_Bits(1) > 0 
%         Re = 1/sqrt(2);
%     else
%         Re =-1/sqrt(2);
%     end
%     if Tx_Bits(2) > 0 
%         Im = 1/sqrt(2);
%     else
%         Im =-1/sqrt(2);
%     end
%     Tx_Symbol = complex(Re, Im);
%     %%% QPSK에 대해 여러분이 코딩을 해야됨!
% end


    if M==2
        figure(100);
        if Tx_Bits == 1
            plot(real(rx_signal), imag(rx_signal), 'bo');
        else
            plot(real(rx_signal), imag(rx_signal), 'rx');
        end
        hold on;
        grid on;
        %legend('Tx bit: 1', 'Tx bit: 0', 'Location','southwest');
        xlabel('Real of Rx Signal', 'Fontsize', 12);
        ylabel('Imag of Rx Signal', 'Fontsize', 12);
        title('Constellation for BPSK', 'Fontsize', 14);
        axis([-3 3 -3 3]);
        %axis([-10 10 -10 10]);
        pause();
    elseif M==4
        figure(400);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   [Tx_Bits(1),Tx_Bits(2)]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   [0,1] | [1,1]
        % --------------------------
        %   [0,0] | [1,0]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Tx_Bits(1) == 1 && Tx_Bits(2) == 1 %%%%%%%%%%%%%%%%%%%%% 1사분면 !!
            plot(real(rx_signal), imag(rx_signal), 'bo');
        elseif Tx_Bits(1) == 0 && Tx_Bits(2) == 1 %%%%%%%%%%%%%%%%% 2사분면 !!
            plot(real(rx_signal), imag(rx_signal), 'rx');
        elseif Tx_Bits(1) == 0 && Tx_Bits(2) == 0 %%%%%%%%%%%%%%%%% 3사분면 !!
            plot(real(rx_signal), imag(rx_signal), 'gs');
        else%if Tx_Bits(1) == 1 && Tx_Bits(2) == 0 %%%%%%%%%%%%%%%% 4사분면 !!
            plot(real(rx_signal), imag(rx_signal), 'k+');
        end
        hold on;
        grid on;
        %legend('Tx bit=[1,1]', 'Tx bit=[0,1]', 'Tx bit=[0,0]', 'Tx bit=[1,0]', 'Location','southwest');
        xlabel('Real of Rx Signal', 'Fontsize', 12);
        ylabel('Imag of Rx Signal', 'Fontsize', 12);
        title('Constellation for QPSK', 'Fontsize', 14);
        axis([-1.5 1.5 -1.5 1.5]);
        %axis([-3 3 -3 3]);
        pause();
    end
    out = 1;
end

