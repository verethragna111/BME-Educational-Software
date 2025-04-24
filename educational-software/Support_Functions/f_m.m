function mp = f_m(m,v_m)
    % Differential Equation of m
    
    mp = alpha_m(v_m).*(1-m) - beta_m(v_m).*m;
    end