function c_hat = get_time_average(params)
    N = params.N;
    GP = params.GP;
    L = params.L;
    r = params.rx_signal;    
    N_OFDM_symbols = length(r)/(GP+N);
    
    c = 0;
    for k = 1+GP:N
        for m = 1:N_OFDM_symbols
            c = c + abs(r(k+(GP+N)*(m-1)))^2;
        end
    end
    c_hat = c/(N_OFDM_symbols*(N-GP)); %%c_hat을 생성

end
