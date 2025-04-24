function rate = beta_h(v)
    rate = zeros(size(v));
    rate =  1.0 ./ (exp((-v+30.0)/10.0) + 1.0);
end 