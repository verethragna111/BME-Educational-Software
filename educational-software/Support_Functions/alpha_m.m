function rate = alpha_m(v)
    rate = zeros(size(v));		% DEFAULT RATE TO ZERO
    rate = 0.1.*(25.-v)./(exp((25.-v)./10)-1);
end 