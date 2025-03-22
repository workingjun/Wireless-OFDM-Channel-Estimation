function [mn_sol,zs_sol,rob,rak_v2] = new_method4(new_maxium,SIM)
    new_maxium(16) = [];%%%%%method2와 method3의  u가 각각 15와 16이기 때문에 method3의 16번째 u를 배제한다.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%minmax-normalization
    mn1 = (new_maxium-min(new_maxium))/(max(new_maxium)-min(new_maxium));
    mn2 = (SIM - min(SIM))/(max(SIM)-min(SIM));
    [~,mn_sol] = max(mn1+mn2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%zscore-normalization
    avg1 = mean(new_maxium);
    sig1 = mean(new_maxium.^2)-mean(new_maxium).^2;
    avg2 = mean(SIM);
    sig2 = mean(SIM.^2)-mean(SIM).^2;
    zs1 = (new_maxium-avg1)/sig1;
    zs2 = (SIM - avg2)/sig2;
    [~,zs_sol] = max(zs1+zs2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%robust-normalization
    rb1 = sort(new_maxium);
    rb2 = sort(SIM);
    rb11 = (new_maxium - median(new_maxium))/(rb1(11)-rb1(4));
    rb21 = (SIM - median(SIM))/(rb2(11)-rb2(4));
    [~,rob] = max(rb11 + rb21);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ranking-sum
    [~,bin]=sort(SIM,'descend');
            
    [~,bin2]=sort(new_maxium,'descend');
    for s = 1:15
        r1(s) = s;
        r2(s) = s;
    end
    for r1 = 1:15
        for r2 = 1:15
            if bin(r1) == bin2(r2)
                r12(r1) = r1 + r2;
            end
        end
    end
    [~,nd] = min(r12);
    rak_v2 = bin(nd);
end