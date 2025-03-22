function [J,SIM,new_sol,e, AVG, sigma] = new_method2(r,GP,N)
        J = zeros(1,16);        %%%% 0으로 이루어진 1x16 행렬 생성
        AVG = zeros(1,16);      %%%% 0으로 이루어진 1x16 행렬 생성
        sigma = zeros(1,16);    %%%% 0으로 이루어진 1x16 행렬 생성
        e = zeros(1,15);        %%%% 0으로 이루어진 1x15 행렬 생성
        %S = ones(1,16);         %%%% 0으로 이루어진 1x15 행렬 생성
        SIM = zeros(1,15);
        N_OFDM_symbols = length(r)/(GP+N);    %ofdm symbol의 갯수를 설정함
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%noise varience estimators를 생성
        for u = 1:GP
            Ttmp = 0;
            for m = 1:N_OFDM_symbols
                for kk = u:GP
                    Ttmp = Ttmp + abs(r((GP+N)*(m-1)+kk)-r((GP+N)*(m-1)+N+kk))^2;
                end
            end
            J(u) = Ttmp/(2*N_OFDM_symbols*(GP-u+1));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for u = 1:GP-1
            e(u) = J(u) - (1-1/(GP-u+1))*J(u+1);
        end

        for u = 1: GP-1
            AVG(u) = J(u) ./ (GP-u+1);             %%%%%평균을 계산
            sigma(u) = (AVG(u).^2)/N_OFDM_symbols; %%%%%분산을 계산
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%가우시안 pdf에서 ML값을 추정한다.
        for u = 1:GP-1
            for mm = u:GP-1
                SIM(u) = SIM(u) - abs(e(mm)-AVG(mm))^2/(2*sigma(mm))-log(2*pi*sigma(mm))/2;
            end
            SIM(u) = SIM(u)/(GP-u);
        end
        [~,new_sol] = max(SIM);
    end

