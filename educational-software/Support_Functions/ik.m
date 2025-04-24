function I_K = ik(v,n,g_K,e_K)

    g_K = 36;
    e_K = -12.26;
    I_K = g_K*n.^4.*(v-e_K);
    end