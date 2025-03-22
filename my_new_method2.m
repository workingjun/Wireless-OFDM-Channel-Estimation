function [J, e_sum_pdf, my_method2_sol] = my_new_method2(params)
    N = params.N;
    GP = params.GP;
    L = params.L;
    r = params.rx_signal;
    N_OFDM_symbols = length(r)/(GP+N);
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mean_e_og = params.J ./ (GP - (1:GP) + 1);
    sigma_e_og = mean_e_og.^2 / N_OFDM_symbols;

    e_logpdf = -0.5 * (abs(params.ABdiffsq - mean_e_og).^2 ./ sigma_e_og) - log(sqrt(2 * pi * sigma_e_og));

    e_sum_pdf = cumsum(e_logpdf, 'reverse');
    e_sum_pdf = e_sum_pdf .* (1 ./ (GP + 1 - (1:GP)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, my_method2_sol] = max(e_sum_pdf);

end
