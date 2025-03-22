%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

params.M=2;% BPSK
BPS = log2(params.M);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bit per Symbols :: �ϳ��� ���� �ɺ��� ��Ʈ �� !!

params.N = 64;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT_size = N;
params.GP = 16;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N+GP => 1-OFDM length !!
params.L = 4;
params.N_OFDM_symbols = 10^2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Symbols :: FFT_size <= N by ��պ� at 201012!
N_bits = (BPS*params.N)*params.N_OFDM_symbols;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Bits per OFDM Symbol == BPS*N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.SNR_dB = 10;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constellation �׸��� Enable !!
% params.SNR_dB = -10:2.5:40;
% params.SNR_dB = -10:5:30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER Simulation !!
params.N_iter = 10^2;%%% ���� ���� ��Ȯ���� ���̷��� ���� Ű��ÿ�! :: ���� ����� 10^5 �̻� ���� !! 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

No_Pilot_symbols = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Number of Pilots OFDM Symbols !! [201119]
if No_Pilot_symbols == 0
    Pilot_CHE_Test = 0;
else
    Pilot_CHE_Test = 1;
end 

%%% SNR (Linear Value, ������) => SNR_dB (dB Value, dB ��)
%%% if SNR = 10, SNR_dB <= 10*log10( SNR ) = 10 (dB). 0.
%%% if SNR =100, SNR_dB <= 10*log10( SNR ) = 20 (dB). 
%%%     .^ :: '.' ��ȣ�� �ִ� �Ͱ� ���� ���� ���� ?
%%% SNR_dB �� ���� Ȥ�� Array ����! => SNR �� ���� Ȥ�� Array ����!  
params.SNR = 10.^(params.SNR_dB/10);%% Define of Linear SNR :: SNR_dB = 10*log10( SNR )

%%% SNR = 1/N_power :: SNR = ��ȣ�������� = ��ȣ_����/����_���� = 1/����_�л�
N_power = 1./params.SNR;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���� ���� �Ŀ�(����) = ���� �л� = sigma^ 2

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

Sim_BER = zeros(1,length(params.SNR_dB));%%% C ���� �迭�� ��� �Ͱ� ����! :: BER �̶�� ������ (1��, length(SNR_dB)��)�� �������·� ��� �ʱⰪ�� '0=zero'�� �Ҵ�! 
Sim_SER = zeros(1,length(params.SNR_dB));

Sim_BER2 = zeros(1,length(params.SNR_dB));%%% C ���� �迭�� ��� �Ͱ� ����! :: BER �̶�� ������ (1��, length(SNR_dB)��)�� �������·� ��� �ʱⰪ�� '0=zero'�� �Ҵ�! 
Sim_SER2 = zeros(1,length(params.SNR_dB));

%%
for n = 1:length(params.SNR_dB)
    %%%%%%%%% ';'��ȣ�� ������ �ش簪�� ȭ�鿡 ���! => ���α׷��� ���ư��� ������ Ȯ���ϱ� �ƹ��� ������ ���� �ʰ� ���� ȭ�鿡 ���� ���!
    params.SNR_dB(n)
    params.SNR_dB_fixed = params.SNR_dB(n);
    params.n = n;

    Bit_Error = 0;%%%%%%%%%%%%%%%%%%% �ʱⰪ!
    Symbol_Error = 0;%%%%%%%%%%%%%%%% �ʱⰪ!

    Bit_Error2 = 0;%%%%%%%%%%%%%%%%%%% �ʱⰪ!
    Symbol_Error2 = 0;%%%%%%%%%%%%%%%% �ʱⰪ!
    
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �۽ű� Part(1) :: ���� ��Ʈ ���� !! 
        %%% randi() :: ������ integer(����) ������ �߻���Ű�� �Լ� !!
        Tx_Bits = randi([0 1],1,N_bits);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randdom Bit Generation ==> �����ϰ����ϴ� ��Ʈ�� !!
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �۽ű� Part(2) :: ���� OFDM �ɺ� ���� !! 
        for tt=1:params.N_OFDM_symbols
            ST = (tt-1)*(BPS*params.N) + 1; 
            To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
            tx_bits_per_OFDM = Tx_Bits(ST:To);%% OFDM �ɺ��� (BPS*N) bits !!
            
            for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N ���� PSK �ɺ� => 1-OFDM symbol !! 
                ST = (ss-1)*BPS + 1; 
                To = (ss-1)*BPS + BPS;
                Signal_F(ss) = Tx_Symbol_gen(params.M, tx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> ������ ���� �ɹ� ���� !!
            end
            ST = (tt-1)*params.N + 1; 
            To = (tt-1)*params.N + params.N; 
            Tx_Symbols(ST:To) = Signal_F;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Symbol_Error ����� ���� Tx_Symbols �� ������!!
            
            Signal_F = fftshift(Signal_F);
            Tx_Signal_T = ifft(Signal_F)*sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Signal_F(���ļ���) => Tx_Signal_T(�ð���) ��ȣ�� ��ȯ !! IDFT == IFFT 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [�������� �߰� !!]

            Tx_Add_GP(1:params.GP) = Tx_Signal_T(params.N-(params.GP-1):params.N);%%%%%%%%%%%%%%%%%%%%% CP ������ GP => ��ü N ���� �ڿ��� GP ���� ���� !!
            Tx_Add_GP(params.GP+1:params.N+params.GP) = Tx_Signal_T;%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tx_Add_GP �� ���� => [GP+N]=[16+64] !! 
            
            ST = (tt-1)*(params.N+params.GP) + 1; 
            To = (tt-1)*(params.N+params.GP) + (params.N+params.GP); 
            tx_signal(ST:To) = Tx_Add_GP;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN ä�� ���
        %%% tx_signal :: ���� OFDM ��ȣ(�ɹ�) !!
        %%% N_AWGN    :: ���� AWGN ���� !! :: randn() => (���0, �л�1)�� ����þ� ���� ������ �߻�!! 
        N_AWGN = complex(randn(1,length(tx_signal)), randn(1,length(tx_signal))).*sqrt(N_power(n)/2);%% Generate AWGN : Additive White Gaussian Noise
        h = complex(randn(1,params.L),randn(1,params.L)) / sqrt(2);
        params.h = h/sqrt(params.L);
        h = params.h;

        H= fft(h,params.N);
        H= fftshift(H);

        ch_out = conv(tx_signal,h); 

        tx_signal = ch_out(1:length(tx_signal));
        params.rx_signal = tx_signal + N_AWGN;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN ä���� ����� ���� ��ȣ ���� !!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ļ� ������ ���̵� ä�� ??????????????????????? (��)

        [params.J,SIM,new_sol] = new_method2(params.rx_signal,params.GP,params.N);%%%%%%(u)th estimated noise power
        
        [new_maxium,L_sol,c_hat,p,rho] = new_method3(params.rx_signal,params.GP,params.N,params.J);%%%%%%%%method3 �Լ��� ������ ���
        % disp('�޸���ƽ�ϰ� ���� ���')
        % disp(p);
        %[mn_sol,zs_sol,rob,rak_v2] = method4_nor(new_maxium,SIM); %%%%%%%%method4�� ranking-sum�� normalization 3����
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        [params.a_k, params.b_k, params.g_k, ...
            params.Asq, params.Bsq, ...
            params.ABdiffsq, params.ABdiffsq_ratio] = get_random_var(params); % random variable ����

        [params.J, params.e_sum_pdf, params.my_method2_sol] = my_new_method2(params);
        [params.rho, params.c_hat, my_method3_sol] = my_new_method3(params);
        
        % filtered_roots = roots_in_range(params.a_k, params.b_k, params.c_hat);
        % disp('�Լ��� ���� ���������� ��')
        % disp(filtered_roots);
        % pause;
        % 
        % Subplot_rxSignal(params) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LKH_Analy(params) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [e_sol1, e_sol2,new_sol2] = method2_upgraded(params);
        
        % figure(1);
        % stem(1:16, params.e_rx_logpdf);
        % ylim([-10, 15]);
        % grid on;
        % pause;

        % Subplot_method2_upgrade_LLR(params, 1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subplot_method2_upgrade(params, false); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if Pilot_CHE_Test == 1 
            hat_H = Pilot_CHE(No_Pilot_symbols, params.N, params.GP, params.rx_signal, Tx_Symbols);
        else
            hat_H = H;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(1) :: ���� ��Ʈ ���� == Rx Bit Detection !! 
        for tt=1:params.N_OFDM_symbols
            ST = (tt-1)*(params.N+params.GP) + params.GP+1; 
            To = (tt-1)*(params.N+params.GP) + params.GP+params.N; 
            Rx_Signal_T = params.rx_signal(ST:To);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ð��࿡�� GP�� ������ OFDM �ɺ�!!
            Rx_Signal_F = fft(Rx_Signal_T)/sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ļ� ������ ��ȯ !!
            Rx_Signal_F = fftshift(Rx_Signal_F);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [�������� �߰� !!]

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (��->��)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ä�ο� ���� �ְ� ����? :: AWGN ä���� �ʿ� ����!!
            Rx_Signal_F_DivH = Rx_Signal_F./hat_H;

            for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N ���� PSK �ɺ� => 1-OFDM symbol !! 
                Bits_tmp = Rx_Bits_gen(params.M, Rx_Signal_F_DivH(ss));
                ST = (ss-1)*BPS + 1; 
                To = (ss-1)*BPS + BPS;
                rx_bits_per_OFDM(ST:To) = Bits_tmp;
            end
            ST = (tt-1)*(BPS*params.N) + 1; 
            To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
            Rx_Bits(ST:To) = rx_bits_per_OFDM;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OFDM �ɺ��� (BPS*N) bits !!

            for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N ���� PSK �ɺ� => 1-OFDM symbol !! 
                ST = (ss-1)*BPS + 1; 
                To = (ss-1)*BPS + BPS;
                Signal_F(ss) = Tx_Symbol_gen(params.M, rx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> ������ ���� �ɹ� ���� !!
            end
            ST = (tt-1)*params.N + 1; 
            To = (tt-1)*params.N + params.N; 
            Rx_Symbols(ST:To) = Signal_F;
            Rx_Symbols_forTxSymbols(ST:To) = Rx_Signal_F;
            H_hat = Rx_Symbols_forTxSymbols(ST:To)./Tx_Symbols(ST:To);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % H_hat = fftshift(H_hat);%%%%%%%%%%%%%%%%%%%%%%%%fft�� �� ���� shift�� ���� ifft�� �ϸ� 1,-1,1,-1�� ���� ����� ��ȯ�� ���� ����(x)           
            h_hat = ifft(H_hat)*sqrt(params.N);%�ð������� ��ȯ, �����л��� 1/N*sqrt(N)^2���� �Ŀ� ����, ä�κл���         
            y_i = abs(h_hat).^2; 

            params.Y_u = ((tt-1)*params.Y_u + y_i)/tt; %���������� ���ؼ� M ���

            L1_sums = cumsum(y_i(1:params.N-1)); %1���� u���� summation -> ���������� ǥ��  
            L2_sums = cumsum(y_i(2:params.N), 'reverse') ./ (params.N-1:-1:1); %u+1���� N���� summation, ���

            params.P(1:params.N-1) = ((tt-1) * params.P(1:params.N-1) + L1_sums) / tt; %���������� ���ؼ� M���
            params.Z(1:params.N-1) = ((tt-1) * params.Z(1:params.N-1) + L2_sums) / tt; %���������� ���ؼ� M���

            y_i_cumsum = cumsum(y_i, 'reverse')./(params.N:-1:1); %�Ųٷ� ���� summation          
            params.JJ = ((tt-1) * params.JJ + y_i_cumsum) / tt; %���������� ���ؼ� M���

        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % [params.p_pdf, params.z_pdf, params.p_Log_pdf, params.z_Log_pdf, params.L_sol_p, params.L_sol_z] = hhat_Taps_Noise(params); %�������� P�� Z�� pdf
        % [params.Y_pdf, params.Y_multi_pdf, params.Y_sum_pdf, params.L_sol_y1, params.L_sol_y2] = hhat_Y(params); %�������� Y_l�� pdf, joint pdf, LLR
        % [params.e_pdf, params.e_multi_pdf, params.e_sum_pdf, params.L_sol_e1, params.L_sol_e2] = hhat_epcilon(params); %�������� e_l�� pdf, joint pdf, LLR
        % [params.Pe_multi_pdf, params.Pe_sum_pdf, params.L_sol_pe1, params.L_sol_pe2] = hhat_pEpcilon(params); %�������� Pe_l�� joint pdf, LLR

        % Subplot_hhat(params);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(2) :: Error Count !!
        
        % Bit_Error = Bit_Error + ( sum(Tx_Bits ~= Rx_Bits) );%%%%%%%%%%%%%%% ��ɾ��� �ǹ�, ������ �����ϱ� �ٶ� !!
        % Symbol_Error = Symbol_Error + ( sum(Tx_Symbols ~= Rx_Symbols) );%%% ��ɾ��� �ǹ�, ������ �����ϱ� �ٶ� !! 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % h_hat_forDiv(1:params.L_sol_pe1) = h_hat(1:params.L_sol_pe1);
        % h_hat_forDiv(params.L_sol_pe1 + 1:params.N)=0;
        % H_hat_forDiv = fft(h_hat_forDiv, params.N);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(3) :: ���� ��Ʈ ���� == Rx Bit Detection !!
        if exist('H_hat_forDiv', 'var')
            for tt=1:params.N_OFDM_symbols
                ST = (tt-1)*(params.N+params.GP) + params.GP+1; 
                To = (tt-1)*(params.N+params.GP) + params.GP+params.N; 
                Rx_Signal_T = params.rx_signal(ST:To);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ð��࿡�� GP�� ������ OFDM �ɺ�!!
                Rx_Signal_F = fft(Rx_Signal_T)/sqrt(params.N);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ļ� ������ ��ȯ !!
                Rx_Signal_F = fftshift(Rx_Signal_F);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fftshfit && *sqrt(N) && /sqrt(N) [�������� �߰� !!]
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (��->��)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ä�ο� ���� �ְ� ����? :: AWGN ä���� �ʿ� ����!!
                Rx_Signal_F_DivH = Rx_Signal_F./H_hat_forDiv;
        
                for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N ���� PSK �ɺ� => 1-OFDM symbol !! 
                    Bits_tmp = Rx_Bits_gen(params.M, Rx_Signal_F_DivH(ss));
                    ST = (ss-1)*BPS + 1; 
                    To = (ss-1)*BPS + BPS;
                    rx_bits_per_OFDM(ST:To) = Bits_tmp;
                end
                ST = (tt-1)*(BPS*params.N) + 1; 
                To = (tt-1)*(BPS*params.N) + (BPS*params.N); 
                Rx_Bits2(ST:To) = rx_bits_per_OFDM;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OFDM �ɺ��� (BPS*N) bits !!
        
                for ss=1:params.N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N ���� PSK �ɺ� => 1-OFDM symbol !! 
                    ST = (ss-1)*BPS + 1; 
                    To = (ss-1)*BPS + BPS;
                    Signal_F(ss) = Tx_Symbol_gen(params.M, rx_bits_per_OFDM(ST:To));%% M-PSK Symbol Generation ==> ������ ���� �ɹ� ���� !!
                end
                ST = (tt-1)*params.N + 1; 
                To = (tt-1)*params.N + params.N; 
                Rx_Symbols2(ST:To) = Signal_F;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���ű� Part(4) :: Error Count !!
    
            Bit_Error2 = Bit_Error2 + ( sum(Tx_Bits ~= Rx_Bits2) );%%%%%%%%%%%%%%% ��ɾ��� �ǹ�, ������ �����ϱ� �ٶ� !!
            Symbol_Error2 = Symbol_Error2 + ( sum(Tx_Symbols ~= Rx_Symbols2) );%%% ��ɾ��� �ǹ�, ������ �����ϱ� �ٶ� !! 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        params = Performance_count(params, new_sol, e_sol1, e_sol2, new_sol2); %%%%%%%%%%%%%%%%%%%%%%%%%
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
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

    % params.Sim_BER(n) = Bit_Error/(N_bits*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)�� BER ��� !!
    % params.Sim_SER(n) = Symbol_Error/(params.N*params.N_OFDM_symbols*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)�� SER ��� !!  
    % 
    % params.Sim_BER2(n) = Bit_Error2/(N_bits*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)�� BER ��� !!
    % params.Sim_SER2(n) = Symbol_Error2/(params.N*params.N_OFDM_symbols*params.N_iter);%%%%%%%%%%%%%%%%%%%%%%% SNR_dB(n)�� SER ��� !! 

    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN ���� �м� ���!
    % Analysis_BER = zeros(1,length(SNR_dB));
    % Analysis_SER = zeros(1,length(SNR_dB));
    % if M == 2 
    %     Analysis_BER = 1/2*erfc( sqrt(SNR) );
    %     Analysis_SER = 1/2*erfc( sqrt(SNR) );
    % elseif M == 4 
    %     Analysis_BER = 1/2*erfc( sqrt(SNR/2) );%% �� SNR/2 �ΰ�?
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SER_app = 2*Q( sqrt(2*NSR)*sin(pi/M) );
    %     Analysis_SER = 1*erfc( sqrt(SNR)*sin(pi/M) );%%%%%%%%%%%%%%%%%%%%%%%%%% �ٻ�ȭ�� ��!!!
    % elseif M == 8 
    %     Analysis_SER = 1*erfc( sqrt(SNR)*sin(pi/M) );%%%%%%%%%%%%%%%%%%%%%%%%%% �ٻ�ȭ�� ��!!!
    %     Analysis_BER = (1/log2(M))*Analysis_SER;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ٻ�ȭ�� ��!!!
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ���� �м� ���!

end

Subplot_performance(params); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotBER_SER(params); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%