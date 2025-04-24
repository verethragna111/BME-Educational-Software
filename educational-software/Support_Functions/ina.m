    
function I_Na = ina(v,m,h,g_Na,e_Na)

    g_Na = 120;
    e_Na = 117.56;
    I_Na = g_Na.*m.^3.*h.*(v-e_Na);
    end