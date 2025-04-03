function Vm = goldman(K_out)

    RTF= 28.0;
    
    % Ion concentrations (mM)
    K_in = 300; 
    Na_in = 20;
    Na_out = 450;
    
    % Permeability ratios
    P_K = 100;  
    P_Na = 100;
    
    w = K_out + (P_Na/P_K)*Na_out;

    y = K_in + (P_Na/P_K)*Na_in;

    Vm = RTF * log(w / y);

end


