function vp = f_v(v,m,n,h,I_stim, e_L)
    % This is the differential equation for the membrane potential in 
    % terms of the four currents and specific membrane capacitance C_m.
    
    % Setting values into variables
    C_m = 1; g_Na = 120;  g_K = 36;  e_K = -12.26;  e_Na = 117.56;
    
    % Differential equation
    vp = -1/C_m.*(ina(v,m,h,g_Na,e_Na) + ik(v,n,g_K,e_K) + ...
    il(v,e_L)+ I_stim);
    end