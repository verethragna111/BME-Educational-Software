function hp = f_h(h, v_m)
    % Differential Equation of h
    
    hp = alpha_h(v_m).*(1-h) - beta_h(v_m).*h;
    end