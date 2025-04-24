function rate = beta_n(v)
    rate = zeros(size(v));
    rate = 0.125*exp(-v/80.0);
end