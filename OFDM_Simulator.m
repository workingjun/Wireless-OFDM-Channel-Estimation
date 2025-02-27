%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

params.M=2;% BPSK
BPS = log2(params.M);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bit per Symbols :: 하나의 전송 심볼당 비트 수 !!

params.N = 64;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT_size = N;
params.GP = 16;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N+GP => 1-OFDM length !!
params.L = 4;
params.N_OFDM_symbols = 10^3;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Symbols :: FFT_size <= N by 고균병 at 201012!
N_bits = (BPS*params.N)*params.N_OFDM_symbols;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Bits per OFDM Symbol == BPS*N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params.SNR_dB = 0:5:30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER Simulation !!
% params.SNR_dB = -10:5:40;
params.SNR_dB = 10;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constellation 그리기 Enable !!
params.N_iter = 10^2;%%% 모의 실험 정확도를 높이려면 수를 키우시오! :: 보고서 제출시 10^5 이상 으로 !! 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

No_Pilot_symbols = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Pilots OFDM Symbols !! [201119]
if No_Pilot_symbols == 0
    Pilot_CHE_Test = 0;
else
    Pilot_CHE_Test = 1;
end 

%%% SNR (Linear Value, 선형값) => SNR_dB (dB Value, dB 값)
%%% if SNR = 10, SNR_dB <= 10*log10( SNR ) = 10 (dB). 0.
%%% if SNR =100, SNR_dB <= 10*log10( SNR ) = 20 (dB). 
%%%     .^ :: '.' 기호가 있는 것과 없는 것의 차이 ?
%%% SNR_dB 는 벡터 혹은 Array 값임! => SNR 도 벡터 혹은 Array 값임!  
params.SNR = 10.^(params.SNR_dB/10);%% Define of Linear SNR :: SNR_dB = 10*log10( SNR )

%%% SNR = 1/N_power :: SNR = 신호대잡음비 = 신호_전력/잡음_전력 = 1/잡음_분산
N_power = 1./params.SNR;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 복소 잡음 파워(정의) = 잡음 분산 = sigma^ 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Signal_F  = zeros(1,params.N);
tx_signal = zeros(1,(params.N+params.GP)*params.N_OFDM_symbols);

tx_bits_per_OFDM = zeros(1,params.N*BPS);
rx_bits_per_OFDM = zeros(1,params.N*BPS);
Tx_Bits = zeros(1,N_bits);
Rx_Bits = zeros(1,N_bits);
Tx_Symbols = zeros(1,params.N*params.N_OFDM_symbols);
Rx_Symbols = zeros(1,params.N*params.N_OFDM_symbols);

params.P = zeros(1, params.N-1);
params.Z = zeros(1, params.N-1);
params.Y_u = zeros(1, params.N);
params.JJ = zeros(1, params.N);

Sim_BER = zeros(1,length(params.SNR_dB));%%% C 에서 배열을 잡는 것과 유사! :: BER 이라는 변수를 (1행, length(SNR_dB)열)의 벡터형태로 잡고 초기값을 '0=zero'로 할당! 
Sim_SER = zeros(1,length(params.SNR_dB));

Sim_BER2 = zeros(1,length(params.SNR_dB));%%% C 에서 배열을 잡는 것과 유사! :: BER 이라는 변수를 (1행, length(SNR_dB)열)의 벡터형태로 잡고 초기값을 '0=zero'로 할당! 
Sim_SER2 = zeros(1,length(params.SNR_dB));

%%
for n = 1:length(params.SNR_dB)
    %%%%%%%%% ';'기호가 없으면 해당값을 화면에 출력! => 프로그램이 돌아가는 과정을 확인하기 아무런 동작을 하지 않고 단지 화면에 값을 출력!
    params.SNR_dB(n)
    params.SNR_dB_fixed = params.SNR_dB(n);
    params.n = n;

    Bit_Error = 0;%%%%%%%%%%%%%%%%%%% 초기값!
    Symbol_Error = 0;%%%%%%%%%%%%%%%% 초기값!

    Bit_Error2 = 0;%%%%%%%%%%%%%%%%%%% 초기값!
    Symbol_Error2 = 0;%%%%%%%%%%%%%%%% 초기값!
    
    params.count1 = 0;
    params.count2 = 0;
    params.count3 = 0;
    params.count4 = 0;

    params.count11 = 0;
    params.count12 = 0;
    params.count13 = 0;
    params.count14 = 0;

    params.count21 = 0;
    params.count22 = 0;
    params.count23 = 0;
    params.count24 = 0;

    params.count31 = 0;
    params.count32 = 0;
    params.count33 = 0;
    params.count34 = 0;

    params.count41 = 0;
    params.count42 = 0;
    params.count43 = 0;
    params.count44 = 0;
    
    for it = 1: params.N_iter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 송신기 Part(1) :: 전송 비트 생성 !! 
        %%% randi() :: 랜덤한 integer(정수) 난수를 발생시키는 함수 !!
        Tx_Bits = randi([0 1],1,N_bits);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randdom Bit Generation ==> 전송하고자하는 비트열 !!
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 송신기 Part(2) :: 전송 OFDM 심볼 생성 !! 
        for tt=1:params.N_OFDM_symbols
            ST = (tt-1)*(BPS*params.N) + 1; 
            To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
            tx_bits_per_OFDM = Tx_Bits(ST:To);%% OFDM 심볼당 (BPS*N) bits !!
            
            for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N 개의 PSK 심볼 => 1-OFDM symbol !! 
                ST = (ss-1)*BPS + 1; 
                To = (ss-1)*BPS + BPS;
                Signal_F(ss) = Tx_Symbol_gen(params.M, tx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> 전송할 변조 심벌 생성 !!
            end
            ST = (tt-1)*params.N + 1; 
            To = (tt-1)*params.N + params.N; 
            Tx_Symbols(ST:To) = Signal_F;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Symbol_Error 계산을 위해 Tx_Symbols 에 저장함!!
            
            Signal_F = fftshift(Signal_F);
            Tx_Signal_T = ifft(Signal_F)*sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal_F(주파수축) => Tx_Signal_T(시간축) 신호로 변환 !! IDFT == IFFT 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [여러분이 추가 !!]

            Tx_Add_GP(1:params.GP) = Tx_Signal_T(params.N-(params.GP-1):params.N);%%%%%%%%%%%%%%%%%%%%% CP 형태의 GP => 전체 N 개중 뒤에서 GP 개를 추출 !!
            Tx_Add_GP(params.GP+1:params.N+params.GP) = Tx_Signal_T;%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tx_Add_GP 의 길이 => [GP+N]=[16+64] !! 
            
            ST = (tt-1)*(params.N+params.GP) + 1; 
            To = (tt-1)*(params.N+params.GP) + (params.N+params.GP); 
            tx_signal(ST:To) = Tx_Add_GP;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN 채널 통과
        %%% tx_signal :: 전송 OFDM 신호(심벌) !!
        %%% N_AWGN    :: 복소 AWGN 잡음 !! :: randn() => (편균0, 분산1)인 가우시안 랜덤 변수를 발생!! 
        N_AWGN = complex(randn(1,length(tx_signal)), randn(1,length(tx_signal))).*sqrt(N_power(n)/2);%% Generate AWGN : Additive White Gaussian Noise
        h = complex(randn(1,params.L),randn(1,params.L)) / sqrt(2);
        params.h = h/sqrt(params.L);
        h = params.h;

        H= fft(h,params.N);
        H= fftshift(H);

        ch_out = conv(tx_signal,h); 

        tx_signal = ch_out(1:length(tx_signal));
        params.rx_signal = tx_signal + N_AWGN;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN 채널을 통과한 수신 신호 생성 !!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 주파수 선택적 페이딩 채널 ??????????????????????? (가)
        
        params.J = get_noise_variance(params);
        [params.a_k, params.b_k, params.g_k, params.Asq, params.Bsq, params.ABplussq, params.ABdiffsq, params.ABdiffsq_ch] = get_random_var(params);
        params.c_hat = get_time_average(params);

        % [~, new2_sol] = new_method2(params);%%%%%%(u)th estimated noise power
        [rho2, new3_sol] = new_method3(params.rx_signal, params.GP, params.N, params.J);%%%%%%%%method3 함수를 가져와 사용
        % [mn_sol, zs_sol, rob, rak_v2] = new_method4(new_maxium,SIM); %%%%%%%%method4로 ranking-sum과 normalization 3가지
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        [params.e_sum_pdf, params.my_method2_sol] = my_new_method2(params);
        [params.rho, my_method3_sol] = my_new_method3(params);

        % Subplot_rxSignal(params) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LKH_Analy(params) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [params.p_rx_logpdf1, params.p_rx_logpdf2, params.p_rx_sum_pdf1, params.p_rx_sum_pdf2, ...
          params.e_rx_logpdf1, params.e_rx_logpdf2, params.e_rx_sum_pdf1, params.e_rx_sum_pdf2, ...
          params.pe_rx_sum_pdf1, params.pe_rx_sum_pdf2, ...
          params.p_sol1, params.p_sol2, params.e_sol1, params.e_sol2, params.pe_sol1, params.pe_sol2] = method2_upgraded(params);

        Subplot_method2_upgrade_LLR(params, 1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subplot_method2_upgrade(params, false); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Pilot_CHE_Test == 1 
            hat_H = Pilot_CHE(No_Pilot_symbols, params.N, params.GP, params.rx_signal, Tx_Symbols);
        else
            hat_H = H;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(1) :: 수신 비트 검출 == Rx Bit Detection !! 
        % for tt=1:params.N_OFDM_symbols
        %     ST = (tt-1)*(params.N+params.GP) + params.GP+1; 
        %     To = (tt-1)*(params.N+params.GP) + params.GP+params.N; 
        %     Rx_Signal_T = params.rx_signal(ST:To);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 시간축에서 GP를 제거한 OFDM 심볼!!
        %     Rx_Signal_F = fft(Rx_Signal_T)/sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 주파수 축으로 변환 !!
        %     Rx_Signal_F = fftshift(Rx_Signal_F);
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [여러분이 추가 !!]
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (가->나)
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 채널에 의한 왜곡 보상? :: AWGN 채널은 필요 없음!!
        %     Rx_Signal_F_DivH = Rx_Signal_F./hat_H;
        % 
        %     for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N 개의 PSK 심볼 => 1-OFDM symbol !! 
        %         Bits_tmp = Rx_Bits_gen(params.M, Rx_Signal_F_DivH(ss));
        %         ST = (ss-1)*BPS + 1; 
        %         To = (ss-1)*BPS + BPS;
        %         rx_bits_per_OFDM(ST:To) = Bits_tmp;
        %     end
        %     ST = (tt-1)*(BPS*params.N) + 1; 
        %     To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
        %     Rx_Bits(ST:To) = rx_bits_per_OFDM;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OFDM 심볼당 (BPS*N) bits !!
        % 
        %     for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N 개의 PSK 심볼 => 1-OFDM symbol !! 
        %         ST = (ss-1)*BPS + 1; 
        %         To = (ss-1)*BPS + BPS;
        %         Signal_F(ss) = Tx_Symbol_gen(params.M, rx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> 전송할 변조 심벌 생성 !!
        %     end
        %     ST = (tt-1)*params.N + 1; 
        %     To = (tt-1)*params.N + params.N; 
        %     Rx_Symbols(ST:To) = Signal_F;
        %     Rx_Symbols_forTxSymbols(ST:To) = Rx_Signal_F;
        %     H_hat = Rx_Symbols_forTxSymbols(ST:To)./Tx_Symbols(ST:To);
        % 
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        %     % H_hat = fftshift(H_hat);%%%%%%%%%%%%%%%%%%%%%%%%fft를 한 이후 shift를 거쳐 ifft를 하면 1,-1,1,-1과 같이 제대로 변환이 되지 않음(x)           
        %     h_hat = ifft(H_hat)*sqrt(params.N);%시간축으로 변환, 잡음분산은 1/N*sqrt(N)^2으로 파워 맞춤, 채널분산은         
        %     y_i = abs(h_hat).^2; 
        % 
        %     params.Y_u = ((tt-1)*params.Y_u + y_i)/tt; %순차적으로 더해서 M 평균
        % 
        %     L1_sums = cumsum(y_i(1:params.N-1)); %1부터 u까지 summation -> 누적합으로 표현  
        %     L2_sums = cumsum(y_i(2:params.N), 'reverse') ./ (params.N-1:-1:1); %u+1부터 N까지 summation, 평균
        % 
        %     params.P(1:params.N-1) = ((tt-1) * params.P(1:params.N-1) + L1_sums) / tt; %순차적으로 더해서 M평균
        %     params.Z(1:params.N-1) = ((tt-1) * params.Z(1:params.N-1) + L2_sums) / tt; %순차적으로 더해서 M평균
        % 
        %     y_i_cumsum = cumsum(y_i, 'reverse')./(params.N:-1:1); %거꾸로 누적 summation          
        %     params.JJ = ((tt-1) * params.JJ + y_i_cumsum) / tt; %순차적으로 더해서 M평균
        % 
        % end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % [params.p_pdf, params.z_pdf, params.p_Log_pdf, params.z_Log_pdf, params.L_sol_p, params.L_sol_z] = hhat_Taps_Noise(params); %랜덤변수 P와 Z의 pdf
        % [params.Y_pdf, params.Y_multi_pdf, params.Y_sum_pdf, params.L_sol_y1, params.L_sol_y2] = hhat_Y(params); %랜덤변수 Y_l의 pdf, joint pdf, LLR
        % [params.e_pdf, params.e_multi_pdf, params.e_sum_pdf, params.L_sol_e1, params.L_sol_e2] = hhat_epcilon(params); %랜덤변수 e_l의 pdf, joint pdf, LLR
        % [params.Pe_multi_pdf, params.Pe_sum_pdf, params.L_sol_pe1, params.L_sol_pe2] = hhat_pEpcilon(params); %랜덤변수 Pe_l의 joint pdf, LLR

        % Subplot_hhat(params);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(2) :: Error Count !!
        
        Bit_Error = Bit_Error + ( sum(Tx_Bits ~= Rx_Bits) );%%%%%%%%%%%%%%% 명령어의 의미, 동작을 이해하기 바람 !!
        Symbol_Error = Symbol_Error + ( sum(Tx_Symbols ~= Rx_Symbols) );%%% 명령어의 의미, 동작을 이해하기 바람 !! 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % h_hat_forDiv(1:params.L_sol_pe1) = h_hat(1:params.L_sol_pe1);
        % h_hat_forDiv(params.L_sol_pe1 + 1:params.N)=0;
        % H_hat_forDiv = fft(h_hat_forDiv, params.N);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(3) :: 수신 비트 검출 == Rx Bit Detection !!
        if exist('H_hat_forDiv', 'var')
            for tt=1:params.N_OFDM_symbols
                ST = (tt-1)*(params.N+params.GP) + params.GP+1; 
                To = (tt-1)*(params.N+params.GP) + params.GP+params.N; 
                Rx_Signal_T = params.rx_signal(ST:To);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 시간축에서 GP를 제거한 OFDM 심볼!!
                Rx_Signal_F = fft(Rx_Signal_T)/sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 주파수 축으로 변환 !!
                Rx_Signal_F = fftshift(Rx_Signal_F);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [여러분이 추가 !!]
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (가->나)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 채널에 의한 왜곡 보상? :: AWGN 채널은 필요 없음!!
                Rx_Signal_F_DivH = Rx_Signal_F./H_hat_forDiv;
        
                for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N 개의 PSK 심볼 => 1-OFDM symbol !! 
                    Bits_tmp = Rx_Bits_gen(params.M, Rx_Signal_F_DivH(ss));
                    ST = (ss-1)*BPS + 1; 
                    To = (ss-1)*BPS + BPS;
                    rx_bits_per_OFDM(ST:To) = Bits_tmp;
                end
                ST = (tt-1)*(BPS*params.N) + 1; 
                To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
                Rx_Bits2(ST:To) = rx_bits_per_OFDM;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OFDM 심볼당 (BPS*N) bits !!
        
                for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N 개의 PSK 심볼 => 1-OFDM symbol !! 
                    ST = (ss-1)*BPS + 1; 
                    To = (ss-1)*BPS + BPS;
                    Signal_F(ss) = Tx_Symbol_gen(params.M, rx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> 전송할 변조 심벌 생성 !!
                end
                ST = (tt-1)*params.N + 1; 
                To = (tt-1)*params.N + params.N; 
                Rx_Symbols2(ST:To) = Signal_F;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 수신기 Part(4) :: Error Count !!
    
            Bit_Error2 = Bit_Error2 + ( sum(Tx_Bits ~= Rx_Bits2) );%%%%%%%%%%%%%%% 명령어의 의미, 동작을 이해하기 바람 !!
            Symbol_Error2 = Symbol_Error2 + ( sum(Tx_Symbols ~= Rx_Symbols2) );%%% 명령어의 의미, 동작을 이해하기 바람 !! 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            params = Performance_count(params, params.my_method2_sol, new3_sol, my_method3_sol); %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

    end
    
    params.M1_1(n) = params.count1/params.N_iter;
    params.M1_2(n) = params.count2/params.N_iter;
    params.M1_3(n) = params.count3/params.N_iter;
    params.M1_4(n) = params.count4/params.N_iter;
    
    params.M2_1(n) = params.count11/params.N_iter;
    params.M2_2(n) = params.count12/params.N_iter;
    params.M2_3(n) = params.count13/params.N_iter;
    params.M2_4(n) = params.count14/params.N_iter;
    
    params.M3_1(n) = params.count21/params.N_iter;
    params.M3_2(n) = params.count22/params.N_iter;
    params.M3_3(n) = params.count23/params.N_iter;
    params.M3_4(n) = params.count24/params.N_iter;
    
    params.M4_1(n) = params.count31/params.N_iter;
    params.M4_2(n) = params.count32/params.N_iter;
    params.M4_3(n) = params.count33/params.N_iter;
    params.M4_4(n) = params.count34/params.N_iter;
    
    params.M5_1(n) = params.count41/params.N_iter;
    params.M5_2(n) = params.count42/params.N_iter;
    params.M5_3(n) = params.count43/params.N_iter;
    params.M5_4(n) = params.count44/params.N_iter;

    params.Sim_BER(n) = Bit_Error/(N_bits*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)의 BER 계산 !!
    params.Sim_SER(n) = Symbol_Error/(params.N*params.N_OFDM_symbols*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)의 SER 계산 !!  

    params.Sim_BER2(n) = Bit_Error2/(N_bits*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)의 BER 계산 !!
    params.Sim_SER2(n) = Symbol_Error2/(params.N*params.N_OFDM_symbols*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)의 SER 계산 !! 
    
end

Subplot_performance(params); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ber_Performance(params); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%