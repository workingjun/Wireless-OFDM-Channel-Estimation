function[] = LLR_Analy(params)
    GP = params.GP;
    L = params.L;
    LKH1 = zeros(1, GP);
    LKH2 = zeros(1, GP);
    LKH3 = zeros(1, GP);
    LKH4 = zeros(1, GP);
    LKH5 = zeros(1, GP);
    Log_term=zeros(1, GP);
    N_OFDM_symbols = 10^2;

    y = params.y_rx;
    e = params.e_rx;
    J = params.J;
    LKH_constent = -log(2 * pi * (J(L+1)^2)/N_OFDM_symbols ) / 2;

    for u=1:GP
        Log_term(u) = log(factorial(GP - u + 1)) / (GP - u + 1);
    end

    for u=1:GP
        LKH_1 = 0;
        LKH_2 = 0;
        for l=u:GP
            LKH_1 = LKH_1 - abs( y(l) - J(L+1) )^2 / (2 * J(L+1)^2/N_OFDM_symbols);
            LKH_2 = LKH_2 - abs( e(l) - (J(L+1) / (GP - l + 1)) )^2 / (2 * J(L+1)^2/( N_OFDM_symbols * (GP - l + 1)^2 ) );
        end
        
        LKH1(u) = LKH_1/(GP - u + 1) + LKH_constent;
        LKH2(u) = LKH_2/(GP - u + 1) + LKH_constent;
        LKH3(u) = LKH2(u) + Log_term(u); 
    end
    
    [~, log_max] = max(LKH5);

    NN = 1:GP;
    
    fig = figure(101);
    fig.Position = [300, 0, 900, 600];
    sgtitle(['SNR: ' num2str(params.SNR_dB_fixed) 'dB']);

    subplot(2, 2, 1);
    stem(NN, LKH1); 
    xlim([0 17]);
    ylim([-5 10])
    title("Y의 LKH");
   
    subplot(2, 2, 2); 
    stem(NN, LKH2);
    xlim([0 17]);
    ylim([-5 10])
    title('Epcilon의 LKH');
    
    subplot(2, 2, 3); 
    stem(NN, LKH3);
    xlim([0 17]);
    ylim([-5 10])
    title(['E의 LKH, log term / Max: '  num2str(log_max)]);

    subplot(2, 2, 4); 
    stem(NN, Log_term);
    xlim([0 17]);
    ylim([-5 10])
    title('Log term');

    pause;
end