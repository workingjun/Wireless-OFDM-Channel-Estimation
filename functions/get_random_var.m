function [a_k, b_k, g_k, Asq, Bsq, ABplussq, ABdiffsq, ABdiffsq_ch] = get_random_var(params)
    GP = params.GP;
    L = params.L;
    N = params.N;
    r = params.rx_signal;
    N_OFDM_symbols = length(r)/(GP+N);

    a_k = zeros(1,GP);
    b_k = zeros(1,GP);
    g_k = zeros(1,GP);
    Asq = zeros(1, GP);
    Bsq = zeros(1, GP);
    ABdiffsq = zeros(1, GP);
    ABplussq = zeros(1, GP);
    ABdiffsq_ch = zeros(1, GP);

    for kk=1:GP
        sumA = 0;
        sumB = 0;
        sumABdiff = 0;
        sumABplus = 0;
        sumABreal = 0;
        a = 0;
        b = 0;
        for m=1:N_OFDM_symbols
            a = a + abs(r((m-1)*(GP+N)+GP-kk+1))^2 + abs(r((GP+N)*m-kk+1))^2;%%%
            b = b + real(r((m-1)*(GP+N)+GP-kk+1)*conj(r((GP+N)*m-kk+1)));%%%%%%%
            sumA = sumA + abs(r((GP+N)*(m-1) + kk))^2;
            sumB = sumB + abs(r((GP+N)*(m-1) + (kk+N)))^2;
            sumABdiff = sumABdiff + abs(r((GP+N)*(m-1) + kk) - r((GP+N)*(m-1) + (kk+N)))^2;
            sumABplus = sumABplus + abs(r((GP+N)*(m-1) + kk) + r((GP+N)*(m-1) + (kk+N)))^2;
            sumABreal = sumABreal + real(r((GP+N)*(m-1) + kk) * conj(r((GP+N)*(m-1) + (kk+N))));
        end
        a_k(kk) = a/N_OFDM_symbols;
        b_k(kk) = b/N_OFDM_symbols; 

        Asq(kk) = sumA/(N_OFDM_symbols);
        Bsq(kk) = sumB/N_OFDM_symbols;
        ABplussq(kk) = sumABplus/(2*N_OFDM_symbols);

        ABdiffsq(kk) = sumABdiff/(2*N_OFDM_symbols*(GP-kk+1));
        ABdiffsq_ch(kk) = sumABdiff/(2*N_OFDM_symbols);
    end
    
    for kk=1:GP+N
        sum_g = 0;
        for m=1:N_OFDM_symbols
            sum_g = sum_g + abs(r((GP+N)*(m-1) + kk))^2;
        end
        g_k(kk) = sum_g/N_OFDM_symbols;
    end
end