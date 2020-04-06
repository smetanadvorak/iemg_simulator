function result = take_scope(signal, center, half_len)
    
if center < half_len+1
    signal = [zeros(half_len + 1 - center, size(signal,2)); signal];
    center = half_len + 1;
end

if center > size(signal,1) - half_len
    signal = [signal; zeros(center - (size(signal,1)-half_len), size(signal,2))];
end

result = signal(center-half_len : center+half_len);

if size(result, 1) ~= 2*half_len+1
    error('Something''s wrong in this function');
end

end

