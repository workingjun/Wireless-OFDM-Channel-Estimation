function [rho, c_hat, my_method3_sol] = my_new_method3(params)
    GP = params.GP;
    N = params.N;
    r = params.rx_signal;
    N_OFDM_symbols = length(r)/(GP+N);
    
    p = zeros(1,GP);
    p_power = zeros(1,GP);
    rho = zeros(GP,GP);
    new_maxium = zeros(1, GP);

    a_k = params.a_k;
    b_k = params.b_k;
    g_k = params.g_k;
    
    c = 0;
    for k = 1+GP:N
        for m = 1:N_OFDM_symbols
            c = c + abs(r(k+(GP+N)*(m-1)))^2;
        end
    end
    c_hat = c/(N_OFDM_symbols*(N-GP)); %%c_hat을 생성

    for kk=1:GP
        for xx=1:-0.001:0
            if method3_power_func(xx, kk, c_hat, a_k, b_k) < 10^-10
                p(kk) = xx;
                break;
            end
        end
    end

    p_power(1) = p(16)*c_hat;
    
    for kk=2:GP
        rho_a = p(GP-kk+1);
        rho_b = p(GP-kk+2);
        if rho_a > rho_b
            p_power(kk) = (rho_a - rho_b)*c_hat;
        else
            p_power(kk) = 0;
        end
    end
    
    for u=1:GP
        P_u = [p_power(1:u), zeros(1, GP-u)];
        P_u = P_u * (c_hat-params.J(u))/sum(p_power(1:u));

        for kk=1:GP
            for ll=1:GP-kk+1
                rho(u, kk) = rho(u, kk) + P_u(ll)/c_hat;
            end
        end
        new_maxium(u) = -sum((a_k-2*rho(u,:).*b_k)./(c_hat*(1-rho(u,:).^2))+log(1-rho(u,:).^2))/(2*(N+GP))-log(pi)/2-((N*log(c_hat)+sum(g_k(GP+1:N)/c_hat))/(2*(N+GP)));%%%%%%%%논문에 나온 더하는 방식을 사용
    end

    [~,my_method3_sol] = max(new_maxium);

end


function result = method3_power_func(x, k, c_hat, a_k, b_k)

    result = c_hat*x^3 - b_k(k)*x^2 + (a_k(k)-c_hat)*x - b_k(k);

end
