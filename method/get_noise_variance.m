function J = get_noise_variance(params)
    N = params.N;
    GP = params.GP;
    L = params.L;
    r = params.rx_signal;    
    N_OFDM_symbols = length(r)/(GP+N);

    J = zeros(1,16); 
    for u = 1:GP
        Ttmp = 0;
        for m = 1:N_OFDM_symbols
            for kk = u:GP
                Ttmp = Ttmp + abs(r((GP+N)*(m-1)+kk)-r((GP+N)*(m-1)+N+kk))^2;
            end
        end
        J(u) = Ttmp/(2*N_OFDM_symbols*(GP-u+1));
    end

end