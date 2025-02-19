function result = Performance_count(params, varargin)
        n = params.n;        
        L = params.L;
        GP = params.GP;
       
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
         
        result.count1 = params.count1;
        result.count2 = params.count2;
        result.count3 = params.count3;
        result.count4 = params.count4;
        
        result.count11 = params.count11;
        result.count12 = params.count12;
        result.count13 = params.count13;
        result.count14 = params.count14;
        
        result.count21 = params.count21;
        result.count22 = params.count22;
        result.count23 = params.count23;
        result.count24 = params.count24;
        
        result.count31 = params.count31;
        result.count32 = params.count32;
        result.count33 = params.count33;
        result.count34 = params.count34;
        
        result.count41 = params.count41;
        result.count42 = params.count42;
        result.count43 = params.count43;
        result.count44 = params.count44;
        
end