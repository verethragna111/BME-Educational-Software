function rate = alpha_h(v)
    rate = zeros(size(v));
    rate = 0.07*exp(-v/20.0);
end
