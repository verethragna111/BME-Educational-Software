function EK = nernst(K_out)
    
    RTF = 28.0; 
    z = 1; % Charge of potassium ion (K+)
    
    % Ion concentrations (mM)
    K_in = 300;
    
    % Nernst equation
    EK = (RTF / z) * log(K_out / K_in);
end


