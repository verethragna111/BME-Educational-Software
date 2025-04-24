function rate =  beta_m(v)
    rate = zeros(size(v));
    rate = 4.0*exp(-v/18.0);
end