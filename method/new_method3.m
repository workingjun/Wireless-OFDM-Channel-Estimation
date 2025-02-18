function [new_maxium, rho, L_sol] = new_method3(params)
    GP = params.GP;
    N = params.N;
    r = params.rx_signal;
    
    p_break = zeros(1,GP);
    rho = zeros(GP,GP);  %%%% 0으로 이루어진 16x16 행렬 생성
    pp = zeros(GP,GP);   %%%% 0으로 이루어진 16x16 행렬 생성
    pn = zeros(1,GP);    %%%% 0으로 이루어진 16x16 행렬 생성  
    power_u = zeros(1,GP);
    new_maxium = zeros(1,GP);
    N_OFDM_symbols = length(r)/(GP+N);
    
    a_k = params.a_k;
    b_k = params.b_k;
    g_k = params.g_k;
    c_hat = params.c_hat;
    
    for k =1:16%%%%%%%%%%%%%%%%%%%%%%%%%%0과 1사이의 rho를 계산한다.
        for p = 1:-0.001:0
            f = c_hat*p^3-b_k(k)*p^2+(a_k(k)-c_hat)*p-b_k(k);%%from (21)
            if f <10^-10
                p_break(k) = p;
                break; 
            end 
        end
    end
    p = p_break;
    
    p_pw = zeros(1,16);           %%%%%power 벡터 생성
    p_pw(1) = p(16)*c_hat;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% from (16)
    for k = 2:16
        if p(GP-k+1) > p(GP-k+2)
            p_pw(k) = (p(GP-k+1)-p(GP-k+2))*c_hat;
        else
            p_pw(k) = 0;
        end
    end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Table(1) in[19]
    for u = 1:16
        for k = 1:u
           power_u(u,k) = p_pw(k);  %%%%%%u paths power selection from(16)
        end
        %%%%%%channel power normalization
        pn(u) = sum(power_u(u,:));%%%%%total power
        for k = 1:u
           pp(u,k) = power_u(u,k) * (c_hat-params.J(u))./pn(u);
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:GP
            for l = 1:GP-k+1
                rho(u,k) = rho(u,k) + pp(u,l)/c_hat;  %%%%%%(u)th rho calculation 
            end
        end
         %%%%%%% L_max를 추정하기 위한 계산
        new_maxium(u) = -sum((a_k-2*rho(u,:).*b_k)./(c_hat*(1-rho(u,:).^2))+log(1-rho(u,:).^2))/(2*(N+GP))-log(pi)/2-(N*log(c_hat)/+(N-GP)/(2*(N+GP)));%%%%%%%%논문에 나온 더하는 방식을 사용

    end
    [~,L_sol] = max(new_maxium);
end