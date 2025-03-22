function [Pe_multi_pdf, Pe_sum_pdf, L_sol_pe1, L_sol_pe2] = hhat_pEpcilon(params)
    N = params.N;    
    Pe_sum_pdf = zeros(1, N);
    Pe_multi_pdf = zeros(1, N);
    
    Pe_multi_pdf(1) = params.e_multi_pdf(1);
    NN_vec = 2:N;
    Pe_multi_pdf(NN_vec) = params.p_pdf(NN_vec-1) .* (params.e_multi_pdf(NN_vec).^(N-NN_vec+1));
    Pe_multi_pdf(NN_vec) = Pe_multi_pdf(NN_vec) .^ (1 ./ (N + 2 - NN_vec));

    Pe_sum_pdf(1) = params.e_sum_pdf(1);
    Pe_sum_pdf(NN_vec) = params.p_Log_pdf(NN_vec-1) + params.e_sum_pdf(NN_vec) .* (N-NN_vec+1);
    Pe_sum_pdf(NN_vec) = Pe_sum_pdf(NN_vec) .* (1 ./ (N + 2 - NN_vec));
    
    [~, L_sol_pe1] = max(Pe_multi_pdf(1:16));
    [~, L_sol_pe2] = max(Pe_sum_pdf(1:16));
    
    L_sol_pe1 = L_sol_pe1 - 1;
    L_sol_pe2 = L_sol_pe2 - 1;
end
