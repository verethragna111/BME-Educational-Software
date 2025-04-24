function rate = alpha_n(v)
    rate = zeros(size(v));
    rate = 0.01.*(10.-v)./(exp((10.-v)/10)-1);
end