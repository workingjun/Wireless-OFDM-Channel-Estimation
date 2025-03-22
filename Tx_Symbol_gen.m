function [Tx_Symbol] = Tx_Symbol_gen(M, Tx_Bits)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
%   [Tx_Bits(1),Tx_Bits(2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [0,1] | [1,1]
% --------------------------
%   [0,0] | [1,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       if M == 2%% Bits(1) 만 사용! 
            %%% S = 2*Bits - 1; %%% Binary Phase Shift Keying Modulation Scheme
            %%% 위 한줄과 아래 4줄이 같은 결과임을 숙지하기 바람! 
            if Tx_Bits(1) == 0
                Tx_Symbol = -1; %%% Binary Phase Shift Keying Modulation Scheme
            else %% if Bits(1) == 1
                Tx_Symbol = complex(1,0);
            end
       elseif M == 4
            if Tx_Bits(1) > 0 
                Re = 1/sqrt(2);
            else
                Re =-1/sqrt(2);
            end
            if Tx_Bits(2) > 0 
                Im = 1/sqrt(2);
            else
                Im =-1/sqrt(2);
            end
            Tx_Symbol = complex(Re, Im);
            %%% QPSK에 대해 여러분이 코딩을 해야됨!
        elseif M == 8
            ind = 4*Tx_Bits(1)+2*Tx_Bits(2)+Tx_Bits(3);
%             Tx_Bits
%             ind
%             pause();
            if ind == 0 
                Tx_Symbol = exp(1j*1/8*pi);
            elseif ind == 1
                Tx_Symbol = exp(1j*3/8*pi);
            elseif ind == 2
                Tx_Symbol = exp(1j*7/8*pi);
            elseif ind == 3
                Tx_Symbol = exp(1j*5/8*pi);
            elseif ind == 4
                Tx_Symbol = exp(-1j*1/8*pi);
            elseif ind == 5
                Tx_Symbol = exp(-1j*3/8*pi);
            elseif ind == 6
                Tx_Symbol = exp(-1j*7/8*pi);
            else %% elseif ind == 7
                Tx_Symbol = exp(-1j*5/8*pi);
            end
            %%% 8PSK에 대해 여러분이 코딩을 해야됨!
        end




end

