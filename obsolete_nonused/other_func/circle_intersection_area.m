function area = circle_intersection_area( c1, c2, R1, R2 )
% Not tested
d = norm(c1 - c2);
if (d < R1 + R2)
    
    a = R1^2;
    b = R2^2;
    
    x = (a - b + d^2) / (2 * d);
    z = x^2;
    y = sqrt(a - z);
    
    if (d < abs(R2 - R1))
        area = pi * min(a, b);
    else
        area = a * asin(y / R1) + b * asin(y / R2) - y * (x + sqrt(z + b - a));
    end
    
else
    area = 0;
end

end

