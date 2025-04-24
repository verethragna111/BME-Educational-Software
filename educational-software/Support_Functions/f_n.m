function np = f_n(n, v_m)
    % Differential Equation of n
    
    np = alpha_n(v_m).*(1-n) - beta_n(v_m).*n;
    end