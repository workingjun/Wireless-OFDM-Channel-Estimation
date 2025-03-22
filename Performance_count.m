function params = Performance_count(params, varargin)       
        L = params.L;
        GP = params.GP;
        
        % L = 4;
        % GP = 16;

        % varargin의 개수를 확인하고, 없는 경우 기본값 할당
        numVarargin = length(varargin);
        
        L_sol1 = []; L_sol2 = []; L_sol3 = []; L_sol4 = []; L_sol5 = [];
        
        if numVarargin >= 1, L_sol1 = varargin{1}; end
        if numVarargin >= 2, L_sol2 = varargin{2}; end
        if numVarargin >= 3, L_sol3 = varargin{3}; end
        if numVarargin >= 4, L_sol4 = varargin{4}; end
        if numVarargin >= 5, L_sol5 = varargin{5}; end

        if L_sol1
            % performance analy
            if L_sol1 == L
                params.count1 = params.count1+1;
                params.count2 = params.count2+1;
             elseif L_sol1 >L && L_sol1 <(L+GP)/2
                 params.count2 = params.count2 + 1;
             elseif L_sol1>=(L+GP)/2
                 params.count3 = params.count3+1;
             elseif L_sol1<L
                 params.count4 = params.count4+1;
             end 
        end

        if L_sol2
             if L_sol2 == L 
                params.count11 = params.count11+1;
                params.count12 = params.count12+1;
             elseif L_sol2 > L  && L_sol2 < (L+GP)/2
                 params.count12 = params.count12 + 1;
             elseif L_sol2 >= (L+GP)/2
                 params.count13 = params.count13+1;
             elseif L_sol2 < L 
                 params.count14 = params.count14+1;
             end
        end
        
         if L_sol3
             if L_sol3 == L
                params.count21 = params.count21+1;
                params.count22 = params.count22+1;
             elseif L_sol3 >L && L_sol3 <(L+GP)/2
                 params.count22 = params.count22 + 1;
             elseif L_sol3>=(L+GP)/2
                 params.count23 = params.count23+1;
             elseif L_sol3<L
                 params.count24 = params.count24+1;
             end 
         end
          
         if L_sol4
             if L_sol4 == L
                params.count31 = params.count31+1;
                params.count32 = params.count32+1;
             elseif L_sol4 >L && L_sol4 <(L+GP)/2
                 params.count32 = params.count32 + 1;
             elseif L_sol4>=(L+GP)/2
                 params.count33 = params.count33+1;
             elseif L_sol4<L
                 params.count34 = params.count34+1;
             end
         end
        
         if L_sol5
             if L_sol5 == L
                params.count41 = params.count41+1;
                params.count42 = params.count42+1;
             elseif L_sol5 >L && L_sol5 <(L+GP)/2
                 params.count42 = params.count42 + 1;
             elseif L_sol5>=(L+GP)/2
                 params.count43 = params.count43+1;
             elseif L_sol5<L
                 params.count44 = params.count44+1;
             end
         end
        
end