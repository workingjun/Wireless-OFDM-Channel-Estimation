function filtered_roots = roots_in_range(a_k, b_k, c_hat)
    % 삼차방정식 ax^3 + bx^2 + cx + d = 0의 근을 구하는 함수
    % 단, 0 <= x <= 1 범위에 있는 근만 반환
    % 계수를 이용하여 근 구하기
    filtered_roots = zeros(1, 16);
    for k=1:16
        all_roots = roots([c_hat, -b_k(k), a_k(k)-c_hat, -b_k(k)]);
        % 실수 근만 선택하고, 0 ~ 1 사이에 있는 근만 필터링
        
        for i=1:length(all_roots)
            if imag(all_roots(i)) == 0 && real(all_roots(i)) > 0 && real(all_roots(i)) < 1
                filtered_roots(k) = all_roots(i);
                break;
            end
        end    
    end
end


