function vec = shift_padding(vec, sh, dim)


vec = circshift(vec, sh, dim);
vec(1:sh) = 0;
vec(end+sh+1:end) = 0;

        

end

